//      ecosystem.cpp
//
//      Copyright 2011 Piotr Bentkowski <bentkowski.piotr@gmail.com>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include "ecosystem.h"
#include "cell.h"
#include "tagging_system.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_log.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>  
#include <string.h>
#include "rngEngine.h"

/**
 * @brief Core method. Constructor - sets an ecosystem, but with no values.
 * 
 */
Ecosystem::Ecosystem() {
}

/**
 * \brief Core method. Sets a new population in given Env conditions and of
 * given size.
 *
 * Generates cells and puts them in a vector
 * (<a href="http://www.cplusplus.com/reference/stl/vector/">STL vector</a>)
 * until reaches the limit of ecosystem capacity defined in defines.h file.
 * Should be called ones in the run (in the very beginning).
 * 
 * @param total_resource_pool - total amount of resource in the whole environment
 * @param Env - current value of the environmental conditions
 * @param minimal_cells_resource - minimal resource below which cell dies
 * @param maximal_resourse_bonus_at_start - maximal bonus a cell cat get at start to be above the minimal limit
 * @param min_genotype_initialized - minimal genotype size at initialisation 
 * @param max_genotype_initialized - maximal genotype at initialisation
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - the width of the env space genes can take
 * @param environment_resolution - resolution of the env space
 * @param one_over_total_env_space_surf - 1/total env surface (const supporting computation
 */
void Ecosystem::setPopulation(double total_resource_pool, double Env,
        double minimal_cells_resource, double maximal_resourse_bonus_at_start,
        int min_genotype_initialized, int max_genotype_initialized,
        double not_super_gene_param, double gene_space_width,
        double environment_resolution, double one_over_total_env_space_surf) {
    EnvFreeResources = total_resource_pool;
    startOfMutationCounter = 0;
    double tempEnv = EnvFreeResources;
    Cell tempCell;
    // allocating a population of cells depending on the capacity of Environment
    while (EnvFreeResources > 0.0) {
        tempEnv = EnvFreeResources;
        PopulationOfCells.push_back(tempCell);
        PopulationOfCells.back().setBrandNewCell(Env, minimal_cells_resource,
                maximal_resourse_bonus_at_start, min_genotype_initialized,
                max_genotype_initialized, not_super_gene_param, gene_space_width,
                environment_resolution, one_over_total_env_space_surf);
        EnvFreeResources = EnvFreeResources - PopulationOfCells.back().CellsResource;
    }
    // preventing the Free Resources from being a negative number
    if (EnvFreeResources < 0.0) {
        EnvFreeResources = tempEnv;
        PopulationOfCells.pop_back();
    }
    numberOfDeaths = 0;
    numberOfBirths = 0;
}

/**
 * \brief Core method. Feeds cells. Has two modes of feeding: (p) proportionally to
 * cell's demand for resource and to resource availability and (s) sequential feeding
 * when cells are aligned in a random queue and fed one by one until depleting resources
 * (not all cells are fed then).
 * 
 * @param mode_of_feeding - a string indicating the mode of feeding
 * @param Env - current value of the environmental conditions
 * @param gene_replication_cost - cost of having one metabolic gene
 * @param core_genes - number of non-reducible genes
 * @param metabolic_costs - constant metabolic cost per cell
 * @param gain_from_the_genotype_ratio - the maximal permitted uptake rate
 */
void Ecosystem::feedCells(std::string mode_of_feeding, double Env,
        double gene_replication_cost, int core_genes, double metabolic_costs,
        double gain_from_the_genotype_ratio) {
    // number of cells in population
    int PopSize = (int) PopulationOfCells.size();
    rngEngine* pRandomNumberGenerator = rngEngine::getInstance();
    double SharedBetweenFactor = 0.0;

    if (mode_of_feeding.compare("p") == 0) {
        // array of the demands of all the cells in this round
        double CellsDemandsInThisRound[PopSize];
        // costs of living - cells lose resource and free resource increase
        for (int k = 0; k < PopSize; k++) {
            // evaluation of internal resources of a cell (just costs)
            PopulationOfCells[k].evaluateInternalResourcesCost(gene_replication_cost,
                    core_genes, metabolic_costs);
            // ageing cells
            PopulationOfCells[k].ageingACell();
            // evaluation of free resources
            EnvFreeResources -= PopulationOfCells[k].recentLoss;
            // zeroing cell's recent gain of resources
            PopulationOfCells[k].zeroRecentGainOfCell();
            // evaluation how much a cell can gain in this round of feeding
            PopulationOfCells[k].evaluateGrowthRateFromGenotypeOnePointEnv(Env);
            // Demand of a cell for resources
            CellsDemandsInThisRound[k] = gain_from_the_genotype_ratio \
                * PopulationOfCells[k].CellsGrossGrowthRateAtTheMoment;
            SharedBetweenFactor += CellsDemandsInThisRound[k];
        }
        SharedBetweenFactor = EnvFreeResources / SharedBetweenFactor;

        // the feeding itself:
        // cells are being fed and Free Resources are being depleted
        for (int j = 0; j < PopSize; j++) {
            //evaluation of the gain of resources of a cell
            PopulationOfCells[j].evaluateInternalResourcesGainProportional
                    (gain_from_the_genotype_ratio, SharedBetweenFactor);
            // evaluation of free resources
            EnvFreeResources -= PopulationOfCells[j].recentGain;
        }
    }
    if (mode_of_feeding.compare("s") == 0) {
        // initialising a queuing array which will decide who will be fed first
        int ShuffleArr[PopSize];
        // filling a queuing array which will decide who will be fed first
        for (int i = 0; i < PopSize; i++) {
            ShuffleArr[i] = i;
        }
        // Shuffling elements of queuing array by randomly exchanging each with one other
        for (int i = 0; i < (PopSize - 1); i++) {
            // Random remaining position
            int r = i + (((int) \
                (10000000.0 * pRandomNumberGenerator->rngMTgetDauble())) % (PopSize - i));
            int temp = ShuffleArr[i];
            ShuffleArr[i] = ShuffleArr[r];
            ShuffleArr[r] = temp;
        }
        // costs of living - cells lose resource and free resource increase
        for (int k = 0; k < PopSize; k++) {
            // evaluation of internal resources of a cell (just costs)
            PopulationOfCells[k].evaluateInternalResourcesCost(gene_replication_cost,
                    core_genes, metabolic_costs);
            // evaluation of free resources
            EnvFreeResources -= PopulationOfCells[k].recentLoss;
            // ageing cells
            PopulationOfCells[k].ageingACell();
            // zeroing cell's recent gain of resources
            PopulationOfCells[k].zeroRecentGainOfCell();
        }
        // the feeding itself:
        // cells are being fed and Free Resources are being depleted
        for (int j = 0; j < PopSize; j++) {
            if (EnvFreeResources >= PopulationOfCells[ShuffleArr[j]].CellsResource) {
                // evaluation how much a cell can gain in this round of feeding
                PopulationOfCells[ShuffleArr[j]].evaluateGrowthRateFromGenotypeOnePointEnv(Env);
                //evaluation of the gain of resources of a cell
                PopulationOfCells[ShuffleArr[j]].evaluateInternalResourcesGainSequential
                        (gain_from_the_genotype_ratio);
                // evaluation of free resources
                EnvFreeResources -= PopulationOfCells[ShuffleArr[j]].recentGain;
            } else {
                break;
            }
        }
    }
    if ((mode_of_feeding.compare("p") != 0) && (mode_of_feeding.compare("s") != 0)) {
        std::cout << "Mode can be only p (proportional) or s (sequential)"\
            << std::endl;
    }
}

/**
 * \brief Core method. Reproduces or/and kills cells.
 *
 * Does three things: (a) kill cells 'by accident' (random death), (b) removes form
 * population cells with not enough resources to live (killing weak ones),
 * (c) reproduces those, who are big enough and mutates its offspring if necessary.
 * Massively uses <a href="http://www.cplusplus.com/reference/stl/vector/">STL
 * vector</a> methods.
 * 
 * @param mode_of_feeding - a string indicating the mode of feeding
 * @param Env - current value of the environmental conditions
 * @param duplication_probability - gene duplication probability
 * @param deletion_probability - gene deletion probability
 * @param mutate_threshold - gene modification probability
 * @param minimum_genotype_size - minimum permitted genotype size
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - the width of the env space genes can take
 * @param environment_resolution - resolution of the env space, 1/total env
 * @param one_over_total_env_space_surf - surface (const supporting computation)
 * @param random_death_rate - probability a cell will randomly dead in this iteration
 * @param minimal_cells_resource - resource below which cell dies
 * @param reproduction_resourse_size - size at which the cell divides 
 * @param disparity_at_reproduction - if in proportional feeding mode this value
 * says how uneven the resource division is during the cell division
 */
void Ecosystem::reproduceOrKillCell(std::string mode_of_feeding, double Env,
        double duplication_probability, double deletion_probability,
        double mutate_threshold, int minimum_genotype_size, double not_super_gene_param,
        double gene_space_width, double environment_resolution,
        double one_over_total_env_space_surf, double random_death_rate,
        double minimal_cells_resource, double reproduction_resourse_size,
        double disparity_at_reproduction) {
    // calling and initialising random number generator
    rngEngine* pRandomNumberGenerator = rngEngine::getInstance();
    numberOfDeaths = 0;
    numberOfBirths = 0;
    int PopSize;
    // kill cells 'by accident' (random death)
    PopSize = (int) PopulationOfCells.size();
    for (int i = PopSize - 1; i >= 0; i--) {
        if ((pRandomNumberGenerator->rngMTgetDauble()) < random_death_rate) {
            EnvFreeResources = EnvFreeResources + PopulationOfCells[i].CellsResource;
            PopulationOfCells.erase(PopulationOfCells.begin() + i);
            // count number of killed
            numberOfDeaths += 1;
        }
    }
    // kill the weak ones
    PopSize = (int) PopulationOfCells.size();
    for (int i = PopSize - 1; i >= 0; i--) {
        if ((PopulationOfCells[i].CellsResource) < minimal_cells_resource) {
            EnvFreeResources = EnvFreeResources + PopulationOfCells[i].CellsResource;
            PopulationOfCells.erase(PopulationOfCells.begin() + i);
            // count number of killed
            numberOfDeaths += 1;
        }
    }
    // let the fittest win! (reproduction and mutation)
    PopSize = (int) PopulationOfCells.size();
    // sequential feeding:
    if (mode_of_feeding.compare("s") == 0) {
        for (int j = 0; j < PopSize; j++) {
            PopulationOfCells[j].setCellReprodFlagFalse();
            if (PopulationOfCells[j].CellsResource > reproduction_resourse_size) {
                // calculating the intake ratio for the reproducing cell
                PopulationOfCells[j].setAgeSinceLastReprod();
                PopulationOfCells[j].setAvarageIntakeTillRepr();
                PopulationOfCells[j].setAvarageTotalGainTillRepr();
                // Dividing the resource between the mother and the offspring
                PopulationOfCells[j].CellsResource = \
                        (PopulationOfCells[j].CellsResource) * 0.5;
                // note the fact this cell has reproduced
                PopulationOfCells[j].setCellReprodFlagTrue();
                PopulationOfCells[j].incrementNumberOfReprod();
                PopulationOfCells[j].setCellsResAfterReprod();
                // create a 'newborn' cell
                PopulationOfCells.push_back(PopulationOfCells[j]);
                // set new cell's params to zero (ones which need to be zeroed)
                PopulationOfCells.back().resetParamsOfNewCell();
                // inheriting mothers tag
                PopulationOfCells.back().MothersTag = PopulationOfCells[j].OwnTag;
                // mutate (if necessary) both cells
                PopulationOfCells[j].mutateCell(Env, duplication_probability,
                        deletion_probability, mutate_threshold, minimum_genotype_size,
                        not_super_gene_param, gene_space_width, environment_resolution,
                        one_over_total_env_space_surf);
                PopulationOfCells.back().mutateCell(Env, duplication_probability,
                        deletion_probability, mutate_threshold, minimum_genotype_size,
                        not_super_gene_param, gene_space_width, environment_resolution,
                        one_over_total_env_space_surf);
                // count number of born to record cell's age at first reproduction
                numberOfBirths += 1;
            }
        }
    }
    // proportional feeding:
    if (mode_of_feeding.compare("p") == 0) {
        double Division_disparity;
        for (int j = 0; j < PopSize; j++) {
            PopulationOfCells[j].setCellReprodFlagFalse();
            if (PopulationOfCells[j].CellsResource > reproduction_resourse_size) {
                PopulationOfCells[j].setAgeSinceLastReprod();
                PopulationOfCells[j].setAvarageIntakeTillRepr();
                PopulationOfCells[j].setAvarageTotalGainTillRepr();
                Division_disparity = 0.5 + disparity_at_reproduction \
                    * pRandomNumberGenerator->rngMTgetDauble();
                // note the fact this cell has reproduced
                PopulationOfCells[j].setCellReprodFlagTrue();
                PopulationOfCells[j].incrementNumberOfReprod();
                // create a 'newborn' cell
                PopulationOfCells.push_back(PopulationOfCells[j]);
                PopulationOfCells.back().CellsResource = (1.0 - Division_disparity) \
                    * (PopulationOfCells[j].CellsResource);
                PopulationOfCells[j].CellsResource = Division_disparity \
                    * (PopulationOfCells[j].CellsResource);
                PopulationOfCells[j].setCellsResAfterReprod();
                // set new cell's params to zero
                PopulationOfCells.back().resetParamsOfNewCell();
                // inheriting mothers tag
                PopulationOfCells.back().MothersTag = PopulationOfCells[j].OwnTag;
                // mutate (if necessary) both cells
                PopulationOfCells[j].mutateCell(Env, duplication_probability,
                        deletion_probability, mutate_threshold, minimum_genotype_size,
                        not_super_gene_param, gene_space_width, environment_resolution,
                        one_over_total_env_space_surf);
                PopulationOfCells.back().mutateCell(Env, duplication_probability,
                        deletion_probability, mutate_threshold, minimum_genotype_size,
                        not_super_gene_param, gene_space_width, environment_resolution,
                        one_over_total_env_space_surf);
                // count number of born
                numberOfBirths += 1;
            }
        }
    }
}

/**
 * \brief Core method. Updates growth rates of all the cells in the ecosystem.
 *
 * Called when environmental conditions change. Updates growth rates of all
 * the cells in the ecosystem basing on theirs genotypes properties. For the
 * sake of model's computational complexity should not be called in regular
 * run (there are better way to update that before feeding).
 * @param Env - current value of the environmental conditions.
 * @return void
 */
void Ecosystem::updateGrowthRatesOfCellsInPopulation(double Env) {
    int PopSize = (int) PopulationOfCells.size();
    for (int i = 0; i < PopSize; i++) {
        PopulationOfCells[i].evaluateGrowthRateFromGenotypeOnePointEnv(Env);
    }
}

/**
 * \brief Core method. Does the horizontal gene transfer.
 *
 * Called at each time step. Iterates through the population and for each cell
 * in a population tosses if it's gonna be a gene donor. Then picks a cell which
 * gonna be a gene acceptor and in the last step it tosses for each gene in the
 * donor genotype if it will get transfered to the acceptors genotype. Uses
 * <a href=http://www.cplusplus.com/reference/clibrary/cmath/floor/>cmath floor()</a>
 * for number rounding.
 * 
 * @param hgt_donor_prob - probability of being a donor
 * @param hgt_gene_trans - of one gene to get transfered
 */
void Ecosystem::horizontalGeneTransfer(double hgt_donor_prob,
    double hgt_gene_trans) {
    int LystSize = PopulationOfCells.size();
    int k;
    bool did_hgt_happened;
    // calling and initialising random number generator
    rngEngine* pRandomNumberGenerator = rngEngine::getInstance();
    Tagging_system* pTagging_system = Tagging_system::getInstance();
    for (int i = 0; i < LystSize; i++) {
        // tossing if a cell will be a donor
        if (pRandomNumberGenerator->rngMTgetDauble() < hgt_donor_prob){
            // picking an acceptor cell
            do {
                k = (int) floor(pRandomNumberGenerator->rngMTgetDauble() \
                    * LystSize);
            } while ( k == i);
            int GenoypeSize = PopulationOfCells[i].CellsGenotype.GenotypeValues.size();
            did_hgt_happened = false;
            for (int j = 0; j < GenoypeSize; j++){
                // tossing if a gene gets transfered 
                if (pRandomNumberGenerator->rngMTgetDauble() < hgt_gene_trans){
                    PopulationOfCells[k].CellsGenotype.GenotypeValues.push_back
                    (PopulationOfCells[i].CellsGenotype.GenotypeValues[j]);
                    // incrementing the counter of accepted genes
                    PopulationOfCells[k].CellsGenotype.TotalHGT_number += 1;
                    did_hgt_happened = true;
                }
                if (did_hgt_happened == true){
                    PopulationOfCells[k].OwnTag = pTagging_system->getTag();
                }
            }
        } 
    }  
}

/**
 * \brief Data collecting method. Prints a number of genes in genotypes of all Cells.
 *
 * Test function. Not called in 'regular' runs.
 */
void Ecosystem::printNumberOfGenesInCells() {
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        PopulationOfCells[i].CellsGenotype.printNumberOfGenes();
    }
}

/**
 * \brief Data collecting method. Returns the mean number of genes in
 * genotypes of all the cells.
 *
 * Uses procedures from <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Mean-and-standard-deviation-and-variance.html">
 * GNU Scientific Library</a> to calculate those means.
 * 
 * @return - mean number of genes (double)
 */
double Ecosystem::returnMeanNumberOfGenesInCells() {
    int LystSize = PopulationOfCells.size();
    double GenotypeLength[LystSize];
    for (int i = 0; i < LystSize; i++) {
        GenotypeLength[i] = PopulationOfCells[i].CellsGenotype.returnNumberOfGenes();
    }
    return gsl_stats_mean(GenotypeLength, 1, LystSize);
}

/**
 * \brief Data collecting method. Returns the standard deviation of the number
 * of genes in genotypes of all the cells.
 *
 * Uses procedures from <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Mean-and-standard-deviation-and-variance.html">
 * GNU Scientific Library</a> to calculate those means.
 * 
 * @return - STD of the number of genes (double)
 */
double Ecosystem::returnStdOfGenesInCells() {
    int LystSize = PopulationOfCells.size();
    double GenotypeLength[LystSize];
    for (int i = 0; i < LystSize; i++) {
        GenotypeLength[i] = PopulationOfCells[i].CellsGenotype.returnNumberOfGenes();
    }
    return gsl_stats_sd(GenotypeLength, 1, LystSize);
}

/**
 * \brief Data collecting method. Returns the mean age of all cells
 * in ecosystem.
 *
 * Uses procedures from <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Mean-and-standard-deviation-and-variance.html">
 * GNU Scientific Library</a> to calculate those means.
 * 
 * @return - mean age of cells (double)
 */
double Ecosystem::returnMeanAgeOfCells() {
    int LystSize = PopulationOfCells.size();
    double AgeOfAllCells[LystSize];
    for (int i = 0; i < LystSize; i++) {
        AgeOfAllCells[i] = PopulationOfCells[i].CellsAge;
    }
    return gsl_stats_mean(AgeOfAllCells, 1, LystSize);
}

/**
 * \brief Data collecting method. Returns the standard deviation of the
 * age of cells.
 *
 * Uses procedures from <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Mean-and-standard-deviation-and-variance.html">
 * GNU Scientific Library</a> to calculate those means.
 * 
 * @return - STD of the age of cells (double)
 */
double Ecosystem::returnSTDAgeOfCells() {
    int LystSize = PopulationOfCells.size();
    double AgeOfAllCells[LystSize];
    for (int i = 0; i < LystSize; i++) {
        AgeOfAllCells[i] = PopulationOfCells[i].CellsAge;
    }
    return gsl_stats_sd(AgeOfAllCells, 1, LystSize);
}

/**
 * \brief Data collecting method. Returns returns the sum of all resources
 * allocated in all the cells.
 *
 * @return - the sum of resources in the cells (double)
 */
double Ecosystem::returnsSummOfCellResources() {
    int LystSize = PopulationOfCells.size();
    double Summ = 0.0;
    for (int i = 0; i < LystSize; i++) {
        Summ = Summ + PopulationOfCells[i].CellsResource;
    }
    return Summ;
}

/**
 * \brief Data collecting method. Writes to file a histogram of sizes of genotypes
 * in given time step.
 *
 * File name is fixed in this method what obviously limits the number of model
 * executions in one directory.
 * @param genome_size_histogram_size - maximal size of genome you expect the simulation
 * will produce (will be used to produce time series of histograms)
 */
void Ecosystem::histGenomeSize(int genome_size_histogram_size) {
    // setting up histogram for collecting data of genome sizes
    gsl_histogram * hist_genome = gsl_histogram_alloc(genome_size_histogram_size);
    double ranges_genome[genome_size_histogram_size + 1];
    double Range_genome = 1.0;
    for (int j = 0; j <= genome_size_histogram_size; j++) {
        ranges_genome[j] = Range_genome;
        Range_genome += 1.0;
    }
    gsl_histogram_set_ranges(hist_genome, ranges_genome, genome_size_histogram_size + 1);
    // felling in histogram with genotype size data
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        gsl_histogram_increment(hist_genome,
                (double) PopulationOfCells[i].CellsGenotype.GenotypeValues.size());
    }
    // writing genome data to a file
    std::ofstream histFile;
    // ooop, watch that line below when changing the file name!
    histFile.open("GenomeSizeData.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < genome_size_histogram_size; k++) {
        histFile << gsl_histogram_get(hist_genome, k) << " ";
    }
    histFile << "\n";
    histFile.close();
    gsl_histogram_free(hist_genome);
}

/**
 * \brief Data collecting method. Returns a histogram of ages of cells.

 * Splits into a histogram cells' age data. Writes that to file. File name is
 * fixed in this method what obviously limits the number of model
 * executions in one direData collecting method. ctory.
 * @param HistSize - the number of bins in the histogram of the age structure
 * @return void
 */
void Ecosystem::histCellsAge(int HistSize) {
    // setting up histogram for collecting data of genome sizes
    //  int HistSize = (int) (MAX_GENOTYPE_INITIALIZED + 0.1 * MAX_GENOTYPE_INITIALIZED) - 1;
    gsl_histogram * hist = gsl_histogram_alloc(HistSize);
    double ranges[HistSize + 1];
    double Range = 0.0;
    for (int j = 0; j <= HistSize; j++) {
        ranges[j] = Range;
        Range += 10.0;
    }
    gsl_histogram_set_ranges(hist, ranges, HistSize + 1);
    // felling in histogram with genotype size data
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        gsl_histogram_increment(hist, (double) PopulationOfCells[i].CellsAge);
    }
    // writing genome data to a file
    std::ofstream histFile;
    // ooop, watch that line below when changing the file name!
    histFile.open("CellsAgeData.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < HistSize; k++) {
        histFile << gsl_histogram_get(hist, k) << " ";
    }
    histFile << std::endl;
    histFile.close();
    gsl_histogram_free(hist);
}

/**
 * \brief Data collecting method. Returns the Shannon index value for population
 * diversity. Cells which have even one gene different will be accounted as members
 * of different clonal strains.
 *
 * Calculates Shannon index in a given time. Uses <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Logarithm-and-Related-Functions.html">
 * GNU Scientific Library logarithm function</a>
 * 
 * @return - Shannon index value (double) for population biodiversity. 
 */
double Ecosystem::returnPopulationShannonIndex() {
    int LystSize = (int) PopulationOfCells.size();
    double popSize = (double) PopulationOfCells.size();
    double abundance = 1.0;
    double Summ = 0.0;
//    double speciesRichness = 0.0;
    // zeroing "if already considered" counter
    for (int i = 0; i < LystSize; i++) {
        PopulationOfCells[i].TestIfCounted = false;
    }
    for (int i = 0; i < LystSize; i++) {
        if (PopulationOfCells[i].TestIfCounted == false) {
            abundance = 1.0;
            for (int j = i + 1; j < LystSize - 1; j++) {
                if (i != j && PopulationOfCells[j].TestIfCounted == false \
                        && PopulationOfCells[i].OwnTag == PopulationOfCells[j].OwnTag) {
                    abundance += 1.0;
                    PopulationOfCells[j].TestIfCounted = true;
                }
            }
            // the proper calculations
            Summ = Summ - ((abundance / popSize) * gsl_sf_log(abundance / popSize));
//            speciesRichness += 1.0;
        }
    }
    return Summ;
    // Term with correction seems not to be the case, as this population is a complete
    // set of ALL individuals, not a sample from a larger population.
    //    return Summ - ((speciesRichness - 1.0) / (2.0 * popSize));
}

/**
 * @brief Data collecting method. Returns the Shannon Index for genes in the
 * population. Each gene, regardless in which cell is located, is treated as an
 * independent entity. Uses <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Logarithm-and-Related-Functions.html">
 * GNU Scientific Library logarithm function</a>
 * 
 * @return Shannon Index value for diversity of genes in population.
 */
double Ecosystem::returnGenesShannonIndex() {
    int LystSize = (int) PopulationOfCells.size();
    int GenesInCell = 0;
    std::vector<Gene> AllTheGenesLyst;
    //creating a mega-genotype containing ALL the genes in a population
    for (int q = 0; q < LystSize; q++) {
        GenesInCell = PopulationOfCells[q].CellsGenotype.GenotypeValues.size();
        for (int p = 0; p < GenesInCell; p++) {
            AllTheGenesLyst.push_back(PopulationOfCells[q].CellsGenotype.
                GenotypeValues[p]);
        }
    }
    int GeneCounter = AllTheGenesLyst.size();
    double NumberOfGenes = (double) GeneCounter;
    bool IfCountedLyst[GeneCounter];
    for (int w = 0; w < GeneCounter; w++) {
        IfCountedLyst[w] = false;
    }
    
    double abundance = 1.0;
    double Summ = 0.0;
    
    // calculating the Shannon Index for genes
    for (int i = 0; i < GeneCounter; i++) {
        if (IfCountedLyst[i] == false) {
            abundance = 1.0;
            for (int j = i + 1; j < GeneCounter - 1; j++) {
                if (i != j && IfCountedLyst[j] == false \
                    && AllTheGenesLyst[i].GeneValue[0] == AllTheGenesLyst[j].GeneValue[0] \
                    && AllTheGenesLyst[i].GeneValue[2] == AllTheGenesLyst[j].GeneValue[2] \
                    && AllTheGenesLyst[i].GeneValue[1] == AllTheGenesLyst[j].GeneValue[1] ){
                    abundance += 1.0;
                    IfCountedLyst[j] = true;
                }
            }
            // the proper calculations
            Summ = Summ - ((abundance / NumberOfGenes) * gsl_sf_log(abundance / NumberOfGenes));
        }
    }
    return Summ;
}

/**
 * \brief Data collecting method. Returns a histogram of fractions of total
 * environmental space taken by genotypes.
 * 
 * @param envelope_of_genotype_hist_size - number of bins in the histogram of the
 * fraction of the total env space taken by the genotype
 * @param envelope_of_genotype_hist_range - width of a bin in this histogram
 */
void Ecosystem::histGenotypeEnvelopeSurface(int envelope_of_genotype_hist_size,
        double envelope_of_genotype_hist_range) {
    // setting up histogram for collecting data of frame sizes
    gsl_histogram * hist = gsl_histogram_alloc(envelope_of_genotype_hist_size);
    double ranges[envelope_of_genotype_hist_size + 1];
    double Range = 0.0;
    for (int j = 0; j <= envelope_of_genotype_hist_size; j++) {
        ranges[j] = Range;
        Range += envelope_of_genotype_hist_range;
    }
    gsl_histogram_set_ranges(hist, ranges, envelope_of_genotype_hist_size + 1);
    // felling in histogram with genotype frame size to env surface ratio data
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        gsl_histogram_increment(hist, PopulationOfCells[i].CellsGenotype.EnvSpaceFraction);
    }
    // writing data to a histogram file
    std::ofstream histFile;
    // ooop, watch that line below when changing the file name!
    histFile.open("FrameSizeData.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < envelope_of_genotype_hist_size; k++) {
        histFile << gsl_histogram_get(hist, k) << " ";
    }
    histFile << std::endl;
    histFile.close();
    gsl_histogram_free(hist);
}

/**
 * \brief Data collecting method. Returns a histogram of maximums of genotypes envelopes.
 *
 * Splits into a histogram maximums of envelopes of genotypes. Writes that to file.
 * 
 * @param maximum_of_genotype_hist_size - number of bins in the genome size histogram
 * @param maximum_of_genotype_hist_range -  width of a bin in this histogram
 */
void Ecosystem::histGenotypeEnvelopeMaximum(int maximum_of_genotype_hist_size,
        double maximum_of_genotype_hist_range) {
    // setting up histogram for collecting data
    gsl_histogram * hist = gsl_histogram_alloc(maximum_of_genotype_hist_size);
    double ranges[maximum_of_genotype_hist_size + 1];
    double Range = 0.0;
    for (int j = 0; j <= maximum_of_genotype_hist_size; j++) {
        ranges[j] = Range;
        Range += maximum_of_genotype_hist_range;
    }
    gsl_histogram_set_ranges(hist, ranges, maximum_of_genotype_hist_size + 1);
    // filling in histogram
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        gsl_histogram_increment(hist, PopulationOfCells[i].CellsGenotype.MaxGrowthRate);
    }
    // writing data to a histogram file
    std::ofstream histFile;
    // ooop, watch that line below when changing the file name!
    histFile.open("FrameMaxData.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < maximum_of_genotype_hist_size; k++) {
        histFile << gsl_histogram_get(hist, k) << " ";
    }
    histFile << std::endl;
    histFile.close();
    gsl_histogram_free(hist);
}

/**
 * \brief Data collecting method. Returns a histogram of resources allocated
 * in the cells.
 *
 * Splits into a histogram resource allocated to all cells in the population.
 * Writes that to file.
 *
 * @param resource_alloc_hist_size - number of bins in the histogram
 * @param resource_alloc_hist_range - width of a bin
 * @param resource_alloc_hist_start - most-left value of the histogram
 */
void Ecosystem::histResourceAllocatedInCells(int resource_alloc_hist_size,
        double resource_alloc_hist_range, double resource_alloc_hist_start) {
    // setting up histogram for collecting data
    gsl_histogram * hist = gsl_histogram_alloc(resource_alloc_hist_size);
    double ranges[resource_alloc_hist_size + 1];
    double Range = resource_alloc_hist_start;
    for (int j = 0; j <= resource_alloc_hist_size; j++) {
        ranges[j] = Range;
        Range += resource_alloc_hist_range;
    }
    gsl_histogram_set_ranges(hist, ranges, resource_alloc_hist_size + 1);
    // filling in histogram
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        gsl_histogram_increment(hist, PopulationOfCells[i].CellsResource);
    }
    // writing data to a histogram file
    std::ofstream histFile;
    // ooop, watch that line below when changing the file name!
    histFile.open("RecourceInCells.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < resource_alloc_hist_size; k++) {
        histFile << gsl_histogram_get(hist, k) << " ";
    }
    histFile << std::endl;
    histFile.close();
    gsl_histogram_free(hist);
}

/**
 * \brief Data collecting method. Returns a histogram of resources allocated
 * in the cells.
 *
 * Splits into a histogram resource allocated to all cells in the population.
 * Writes that to file.
 *
 * @param resource_uptake_hist_size - number of bins in the histogram
 * @param resource_uptake_hist_range - width of a bin
 */
void Ecosystem::histResourceUptakenByCells(int resource_uptake_hist_size,
        double resource_uptake_hist_range) {
    // setting up histogram for collecting data
    gsl_histogram * hist = gsl_histogram_alloc(resource_uptake_hist_size);
    double ranges[resource_uptake_hist_size + 1];
    double Range = 0.0;
    for (int j = 0; j <= resource_uptake_hist_size; j++) {
        ranges[j] = Range;
        Range += resource_uptake_hist_range;
    }
    gsl_histogram_set_ranges(hist, ranges, resource_uptake_hist_size + 1);
    // filling in histogram
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        gsl_histogram_increment(hist, PopulationOfCells[i].recentGain);
    }
    // writing data to a histogram file
    std::ofstream histFile;
    // ooop, watch that line below when changing the file name!
    histFile.open("RecourceUptakenCells.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < resource_uptake_hist_size; k++) {
        histFile << gsl_histogram_get(hist, k) << " ";
    }
    histFile << std::endl;
    histFile.close();
    gsl_histogram_free(hist);
}

/**
 * \brief Data collecting method. Splits into a histograms a series of data of
 * all the cells which reproduced.
 *
 * Makes five histograms only for cells which reproduced: the average intake
 * per unit of time of cells; the total gain per cell; cells' age at reproduction;
 * number of reproductions of cells; gene number in cells' genotypes. Writes that
 * to files.
 * 
 * @param at_repr_intake_hist_size - number of bins in the histogram of average intake per cell
 * @param at_repr_intake_bin_range - width of bins in it
 * @param AAR_hist_size - number of bins in age at reproduction histogram 
 * @param AAR_hist_range - width of a bin of age at reproduction histogram 
 * @param Reproduce_hist_size - number of bins in histogram of how many times a cell reproduced
 * @param Reproduce_hist_range - width of a bin of histogram of how many times a cell reproduced
 * @param genome_size_histogram_size - number of bins in genome size histogram
 */
void Ecosystem::histsAgeAndResourceAtReprod(int at_repr_intake_hist_size,
        double at_repr_intake_bin_range, int AAR_hist_size, double AAR_hist_range,
        int Reproduce_hist_size, double Reproduce_hist_range,
        int genome_size_histogram_size) {
    // setting up histogram for collecting data
    gsl_histogram * hist_intake = gsl_histogram_alloc(at_repr_intake_hist_size);
    gsl_histogram * hist_tot_gain = gsl_histogram_alloc(at_repr_intake_hist_size);
    gsl_histogram * hist_age_repr = gsl_histogram_alloc(AAR_hist_size);
    gsl_histogram * hist_num_repr = gsl_histogram_alloc(Reproduce_hist_size);
    gsl_histogram * hist_genome = gsl_histogram_alloc(genome_size_histogram_size);
    double ranges[at_repr_intake_hist_size + 1];
    double ranges_repr[AAR_hist_size + 1];
    double ranges_num_repr[Reproduce_hist_size + 1];
    double ranges_genome[genome_size_histogram_size + 1];
    double Range = 0.0;
    double Range_AAR = 0.0;
    double Range_num_repr = 0.0;
    double Range_genome = 1.0;
    for (int j = 0; j <= at_repr_intake_hist_size; j++) {
        ranges[j] = Range;
        Range += at_repr_intake_bin_range;
    }
    for (int j = 0; j <= AAR_hist_size; j++) {
        ranges_repr[j] = Range_AAR;
        Range_AAR += AAR_hist_range;
    }
    for (int j = 0; j <= Reproduce_hist_size; j++) {
        ranges_num_repr[j] = Range_num_repr;
        Range_num_repr += Reproduce_hist_range;
    }
    for (int j = 0; j <= genome_size_histogram_size; j++) {
        ranges_genome[j] = Range_genome;
        Range_genome += 1.0;
    }
    gsl_histogram_set_ranges(hist_intake, ranges, at_repr_intake_hist_size + 1);
    gsl_histogram_set_ranges(hist_tot_gain, ranges, at_repr_intake_hist_size + 1);
    gsl_histogram_set_ranges(hist_age_repr, ranges_repr, AAR_hist_size + 1);
    gsl_histogram_set_ranges(hist_num_repr, ranges_num_repr, Reproduce_hist_size + 1);
    gsl_histogram_set_ranges(hist_genome, ranges_genome, genome_size_histogram_size + 1);
    // filling in histogram
    int LystSize = PopulationOfCells.size();
    for (int i = 0; i < LystSize; i++) {
        if (PopulationOfCells[i].DidReproduced == true) {
            gsl_histogram_increment(hist_intake,
                    PopulationOfCells[i].AvarageIntakeTillRepr);
            gsl_histogram_increment(hist_tot_gain,
                    PopulationOfCells[i].AvarageTotalGainTillRepr);
            gsl_histogram_increment(hist_age_repr,
                    PopulationOfCells[i].TimeSinceLastReprod);
            gsl_histogram_increment(hist_num_repr,
                    PopulationOfCells[i].NumberOfReprod);
            gsl_histogram_increment(hist_genome,
                    (double) PopulationOfCells[i].CellsGenotype.GenotypeValues.size());
        }
    }
    // writing data to a histogram file
    std::ofstream histFile_intake;
    std::ofstream histFile_tot_gain;
    std::ofstream histFile_repr;
    std::ofstream histFile_num_repr;
    std::ofstream histFile_gene_num;
    // ooop, watch that 4 line below when changing the file name!
    histFile_intake.open("AvarIntakeAtRepr.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    histFile_tot_gain.open("AvarGainAtRepr.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    histFile_repr.open("AgeTillReproduction.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    histFile_num_repr.open("NumberOfReproductions.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    histFile_gene_num.open("GeneNumberOfRepr.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < at_repr_intake_hist_size; k++) {
        histFile_intake << gsl_histogram_get(hist_intake, k) << " ";
        histFile_tot_gain << gsl_histogram_get(hist_tot_gain, k) << " ";
    }
    for (int l = 0; l < AAR_hist_size; l++) {
        histFile_repr << gsl_histogram_get(hist_age_repr, l) << " ";
    }
    for (int m = 0; m < Reproduce_hist_size; m++) {
        histFile_num_repr << gsl_histogram_get(hist_num_repr, m) << " ";
    }
    for (int k = 0; k < genome_size_histogram_size; k++) {
        histFile_gene_num << gsl_histogram_get(hist_genome, k) << " ";
    }
    histFile_intake << std::endl;
    histFile_tot_gain << std::endl;
    histFile_repr << std::endl;
    histFile_num_repr << std::endl;
    histFile_gene_num << std::endl;

    histFile_intake.close();
    histFile_tot_gain.close();
    histFile_repr.close();
    histFile_num_repr.close();
    histFile_gene_num.close();

    gsl_histogram_free(hist_intake);
    gsl_histogram_free(hist_tot_gain);
    gsl_histogram_free(hist_age_repr);
    gsl_histogram_free(hist_num_repr);
    gsl_histogram_free(hist_genome);
}

/**
 * \brief Data collecting method. Returns the mean size of the surface of the
 * environment taken by genotypes envelopes.
 *
 * Uses procedures from
 * <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Mean-and-standard-deviation-and-variance.html">
 * GNU Scientific Library</a> to calculate those means.
 * 
 * @return - mean envelope to environment ratio
 */
double Ecosystem::returnMeanEnvelopeToEnvRatio() {
    int LystSize = PopulationOfCells.size();
    double EnvelopeToEnvRatio[LystSize];
    for (int i = 0; i < LystSize; i++) {
        EnvelopeToEnvRatio[i] = PopulationOfCells[i].CellsGenotype.EnvSpaceFraction;
    }
    return gsl_stats_mean(EnvelopeToEnvRatio, 1, LystSize);
}

/**
 * \brief Data collecting method. Returns the standard deviation of the surface
 * of the environment covered by genotypes envelopes.
 *
 * Uses procedures from <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Mean-and-standard-deviation-and-variance.html">
 * GNU Scientific Library</a> to calculate those means.
 * 
 * @return - STD of envelope to environment ratio
 */
double Ecosystem::returnSTDEnvelopeToEnvRatio() {
    int LystSize = PopulationOfCells.size();
    double EnvelopeToEnvRatio[LystSize];
    for (int i = 0; i < LystSize; i++) {
        EnvelopeToEnvRatio[i] = PopulationOfCells[i].CellsGenotype.EnvSpaceFraction;
    }
    return gsl_stats_sd(EnvelopeToEnvRatio, 1, LystSize);
}

/**
 * @brief Data collecting method. Saving general system's properties to the file.
 *
 * @param k - number of iteration (time stamp)
 * @param env - current value of the environmental conditions
 */
void Ecosystem::writeGeneralDataToFile(int k, double env) {
    std::ofstream GeneralDataFile;
    GeneralDataFile.open("GeneralData.dat", std::ios::out | std::ios::ate | std::ios::app);
    GeneralDataFile << k << " " << env << " " << returnMeanNumberOfGenesInCells() << " " \
        << returnStdOfGenesInCells() << " " << PopulationOfCells.size() << " " \
        << returnsSummOfCellResources() << " " << EnvFreeResources << " " \
        << returnMeanAgeOfCells() << " " << returnSTDAgeOfCells() << " " \
        << returnGenesShannonIndex() << " " << numberOfBirths << " " \
        << numberOfDeaths << " " << returnMeanEnvelopeToEnvRatio() << " " \
        << returnSTDEnvelopeToEnvRatio() << std::endl;
    GeneralDataFile.close();
}

/**
 * @brief Data collecting method. Calculating the speed of evolution.
 *
 * Iterates through the population writing to a file cells' Own Tag and the total
 * numbers of all three kinds of mutation (point mutation, gene duplication gene
 * deletion) which lead to the current genetic lineages. Time stamps must be
 * given by the programmer.
 * 
 * @param current_time_stamp - a time stamp representing the current iteration number
 */
void Ecosystem::cumulativeMutationNumbers(int current_time_stamp) {
    int PopSize = (int) PopulationOfCells.size();
    std::ofstream MutationRecords;
    MutationRecords.open("MutationRecords.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    if(PopSize > 0){
        for (int i = 0; i < PopSize; i++) {
            MutationRecords << PopulationOfCells[i].OwnTag << " " << \
                    PopulationOfCells[i].CellsGenotype.TotalPoinMutationNumber << \
                    " " << PopulationOfCells[i].CellsGenotype.TotalDuplicationNumber << \
                    " " << PopulationOfCells[i].CellsGenotype.TotalDeletionNumber << \
                    " " << PopulationOfCells[i].CellsGenotype.TotalHGT_number << \
                    std::endl;
        }

    } else {
        MutationRecords << "-1 All is dead!" << std::endl;
    }
    MutationRecords << "#start_of_the_mutation_counter: " << startOfMutationCounter << \
        std::endl;
    MutationRecords << "#time_at_data_collecting: " << current_time_stamp << std::endl;
    MutationRecords.close();
}

/**
 * @brief Data collecting method. Clears all three record of all the cells
 * regarding the number of mutations.
 * 
 * Iterates through the population setting all the mutation counters of all three
 * mutation types to zero. Can be called at any time. 
 * @param current_iteration - number of the current iteration.
 */
void Ecosystem::clearMutationCounters(int current_iteration) {
    int PopSize = (int) PopulationOfCells.size();
    if(PopSize > 0){
        for (int i = 0; i < PopSize; i++) {
            PopulationOfCells[i].CellsGenotype.TotalPoinMutationNumber = 0;
            PopulationOfCells[i].CellsGenotype.TotalDuplicationNumber = 0;
            PopulationOfCells[i].CellsGenotype.TotalDeletionNumber = 0;
            PopulationOfCells[i].CellsGenotype.TotalHGT_number = 0;
        }
        startOfMutationCounter = current_iteration;
    }
}

/**
 * @brief Data collecting method. Writes to file all the genotypes for the current
 * time step.
 * 
 * Problem with this method is that the produced data file can go beyond the
 * limit of a line length of a text file. To change that go to file gene.cpp
 * and in function geneToString() change std::setprecision(n) to a lower number.
 * @param ENV - current value of the environmental conditions
 */
void Ecosystem::writeGenotypesToFile(double ENV){
    int LystSize = (int) PopulationOfCells.size();
    double abundance = 1.0;
    // zeroing "if already considered" counter
    for (int i = 0; i < LystSize; i++) {
        PopulationOfCells[i].TestIfCounted = false;
    }
    // writing genome data to a file
    std::ofstream theFile;
    // ooop, watch that line below when changing the file name!
    theFile.open("GenomesOfPopulation.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int i = 0; i < LystSize; i++) {
        if (PopulationOfCells[i].TestIfCounted == false) {
            abundance = 1.0;
            for (int j = i + 1; j < LystSize - 1; j++) {
                if (i != j && PopulationOfCells[j].TestIfCounted == false \
                        && PopulationOfCells[i].OwnTag == PopulationOfCells[j].OwnTag) {
                    abundance += 1.0;
                    PopulationOfCells[j].TestIfCounted = true;
                }
            }
            theFile << abundance << " " \
                    << PopulationOfCells[i].MothersTag << " " << PopulationOfCells[i].OwnTag \
                    << " " << PopulationOfCells[i].CellsGenotype.genotypeToString() \
                    << std::endl;
        }
    }
    theFile.close();
}

/**
 * @brief Data collecting method. Generates 3 files containing data about shape
 * of the population's genotypes.
 * 
 * Generates three data files: Histograms with frequency of \f$U_{i}(x)\f$
 * (see thesis for details); shape of the genotype averaged over the whole
 * population; STD of the genotype averaged over the whole population.
 * 
 * @param environment_resolution - resolution of the env space
 */
void Ecosystem::averageGenomeShape(double environment_resolution){
    int xSpaceSize = (int) PopulationOfCells[0].CellsGenotype.GenomeLandscape.size();
    int PopSize = (int) PopulationOfCells.size();
    double XSpaceMeanArr[xSpaceSize];
    double XSpaceSTDArr[xSpaceSize];
    double CellValueArr[PopSize];
    double theGLvalue;
    
    int GL_HistSize = (int) (1.0 / environment_resolution); //  czy tu nie powinno być 2.0 / envi.... ?
    gsl_histogram * hist = gsl_histogram_alloc(GL_HistSize);
    double ranges[GL_HistSize + 1];
    double GLrange = 0.0;
    for (int j = 0; j <= GL_HistSize; j++) {
        ranges[j] = GLrange;
        GLrange += environment_resolution;
    }
    gsl_histogram_set_ranges(hist, ranges, GL_HistSize + 1);
    
    for (int i = 0; i < xSpaceSize; i++){
        for (int j = 0; j < PopSize; j++){
            theGLvalue = PopulationOfCells[j].CellsGenotype.GenomeLandscape[i];
            CellValueArr[j] = theGLvalue;
            gsl_histogram_increment(hist, theGLvalue);
        }
        XSpaceMeanArr[i] = gsl_stats_mean(CellValueArr, 1, xSpaceSize);
        XSpaceSTDArr[i] = gsl_stats_sd(CellValueArr, 1, xSpaceSize);
    }
    std::ofstream histGenotype;
    histGenotype.open("AvaregeGenotype.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int k = 0; k < GL_HistSize; k++) {
        histGenotype << gsl_histogram_get(hist, k) << " ";
    }
    histGenotype << "\n";
    histGenotype.close();
    gsl_histogram_free(hist);
    
    std::ofstream meanGenotEnvel;
    std::ofstream stdGenotEnvel;
    meanGenotEnvel.open("GenotypeEnvelMean.dat", std::ios::out | std::ios::ate | std::ios::app);
    stdGenotEnvel.open("GenotypeEnvelSTD.dat", std::ios::out | std::ios::ate | std::ios::app);
    for (int m = 0; m < xSpaceSize; m++){
        meanGenotEnvel << std::setprecision(4) << XSpaceMeanArr[m] << " ";
        stdGenotEnvel << std::setprecision(4) << XSpaceSTDArr[m] << " ";
    }
    meanGenotEnvel << "\n";
    meanGenotEnvel.close();
    stdGenotEnvel << "\n";
    stdGenotEnvel.close();
}