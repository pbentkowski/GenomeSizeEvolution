//      cell.cpp
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
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include "cell.h"
#include "tagging_system.h"
#include "rngEngine.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

/**
 * @brief Core method. Constructor - sets a cell, but with no values.
 *
 * Just a constructor.
 */
Cell::Cell() {
}

/**
 * @brief Core method. Sets a brand new cell with initial (random) values.
 *
 * Core method. Sets a brand new cell, randomly assigning its model values and
 * calculates growth rate cell has in given environmental conditions. Also sets
 * initial values of parameters used for data analysis.
 *
 * @param Env - current value of the environmental conditions
 * @param minimal_cells_resource - resource below which cell dies
 * @param maximal_resourse_bonus_at_start - maximal bonus a cell cat get at start to be above the minimal limit
 * @param min_genotype_initialized - minimal genotype size at initialisation 
 * @param max_genotype_initialized - maximal genotype at initialisation
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - the width of the env space genes can take
 * @param environment_resolution - resolution of the env space
 * @param one_over_total_env_space_surf - 1/total env surface (const supporting computation)
 *
 */
void Cell::setBrandNewCell(double Env, double minimal_cells_resource,
        double maximal_resourse_bonus_at_start, int min_genotype_initialized,
        int max_genotype_initialized, double not_super_gene_param,
        double gene_space_width, double environment_resolution,
        double one_over_total_env_space_surf) {
    // calling and initialising random number generator
    rngEngine* p_rngEngine = rngEngine::getInstance();

    Tagging_system* pTagging_system = Tagging_system::getInstance();
    // starting with some resources
    CellsResource = minimal_cells_resource + (maximal_resourse_bonus_at_start \
        * p_rngEngine->rngMTgetDauble());
    recentGain = 0.0;
    // initialising a randomly assigned genotype
    CellsGenotype.setGenotype(min_genotype_initialized, max_genotype_initialized,
            not_super_gene_param, gene_space_width, environment_resolution);
    // calculating the Gross Growth Rate for current conditions
    evaluateGrowthRateFromGenotypeOnePointEnv(Env);
    // a mutation sensor (flag-like)
    TestIfMutated = false;
    // age of a new born cell
    CellsAge = 0.0;
    // just a counter for other functions
    TestIfCounted = false;
    // assigning tags
    OwnTag = pTagging_system->getTag();
    // assigning tags
    MothersTag = OwnTag;
    // zeroing the age at the first reproduction
    AgeAtLastReprod = 0;
    // zeroing the number of reproductions
    DidReproduced = false;
    NumberOfReprod = 0;
    // zeroing the average intake of resource for a cell on its way to reproduce
    AvarageIntakeTillRepr = 0.0;
    AvarageTotalGainTillRepr = 0.0;
    TotalGainSinceRepr = 0.0;
    TotalLoss = 0.0;
    CellsResourceAfterReprod = CellsResource;
    CellsGenotype.evaluateAdaptiveLandcape(environment_resolution,
            one_over_total_env_space_surf);
}

/**
 * @brief Core method. Does mutations do cell's genotype.
 *
 * Mutates some genes using Gene class methods. If cell got mutated calculates
 * new growth rate based on new genotype and assigns new tag for freshly
 * established clone.
 * @param Env -current value of the environmental conditions
 * @param duplication_probability - gene duplication probability
 * @param deletion_probability -gene deletion probability
 * @param mutate_threshold - gene modification probability
 * @param minimum_genotype_size - minimum permitted genotype size
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - he width of the env space genes can take
 * @param environment_resolution - resolution of the env space
 * @param one_over_total_env_space_surf -1/total env surface (const supporting computation)
 *
 */
void Cell::mutateCell(double Env, double duplication_probability,
        double deletion_probability, double mutate_threshold, int minimum_genotype_size,
        double not_super_gene_param, double gene_space_width, double environment_resolution,
        double one_over_total_env_space_surf) {
    // calling and initialising random number generator
    rngEngine* p_rngEngine = rngEngine::getInstance();
    Tagging_system* pTagging_system = Tagging_system::getInstance();
    // this is Cell's class data field, if not mutated will be false
    TestIfMutated = false;
    // calculating what's the probability of any duplication and if positive
    // entering duplication method
    if (p_rngEngine->rngMTgetDauble() < (1.0 - gsl_pow_int(1.0 - duplication_probability,
            (int) CellsGenotype.GenotypeValues.size()))) {
        CellsGenotype.duplicateGene(duplication_probability);
        TestIfMutated = true;
    }
    // calculating what's the probability of any deletion and if positive entering
    // deletion method
    if (p_rngEngine->rngMTgetDauble() < (1.0 - gsl_pow_int(1.0 - deletion_probability,
            (int) CellsGenotype.GenotypeValues.size()))) {
        CellsGenotype.deleteGenes(deletion_probability, minimum_genotype_size);
        TestIfMutated = true;
    }
    // calculating what's the probability of ANY single-gene mutation appearing and if
    // positive, then entering mutation method
    if (p_rngEngine->rngMTgetDauble() > gsl_pow_int(1.0 - mutate_threshold,
            (int) CellsGenotype.GenotypeValues.size())) {
        CellsGenotype.mutateGenes(mutate_threshold, not_super_gene_param,
                gene_space_width);
        TestIfMutated = true;
    }
    if (TestIfMutated == true) {
        //evaluateGrowthRateFromGenotypeOnePointEnv(Env)
        // wasn't necessary to call it here
        OwnTag = pTagging_system->getTag();
    }
    CellsGenotype.evaluateAdaptiveLandcape(environment_resolution,
            one_over_total_env_space_surf);
}

/**
 * @brief Core method. Evaluates cell's growth rate for one given value of
 * environment space.
 *
 * It evaluates growth rate of the cell for one value environment space. Uses
 * <a href="http://www.gnu.org/software/gsl/manual/html_node/Vectors.html">GNU
 * Scientific Library vectors</a> to compute this value.
 * @param Env - current value of the environmental conditions.
 */
void Cell::evaluateGrowthRateFromGenotypeOnePointEnv(double Env) {
    // measuring genome size
    int GenomeSize = (int) CellsGenotype.GenotypeValues.size();
    gsl_vector * TempVector = gsl_vector_alloc(GenomeSize);
    for (int j = 0; j < GenomeSize; j++) {
        // calculating growth rate for a given gene for a current Environment value
        gsl_vector_set(TempVector, j, CellsGenotype.calculateGrowthFromOneGene(j, Env));
    }
    // picking the maximal value for given Env conditions
    CellsGrossGrowthRateAtTheMoment = gsl_vector_max(TempVector);
    gsl_vector_free(TempVector);
}

/**
 * @brief Core method. Evaluates cell's internal resource gain for sequential feeding.
 *
 * Evaluates how much resources cell has gained based on the cell's genotype
 * quality.
 * @param gain_from_the_genotype_ratio - the maximal permitted uptake rate.
 *
 */
void Cell::evaluateInternalResourcesGainSequential(double gain_from_the_genotype_ratio) {
    double recentGainTemp;
    // gain for a cell in one iteration
    recentGainTemp = gain_from_the_genotype_ratio * CellsGrossGrowthRateAtTheMoment;
    recentGain = recentGainTemp;
    CellsResource += recentGainTemp;
    TotalGainSinceRepr += recentGainTemp;
}

/**
 * @brief Core method. Evaluates cell's internal resource gain for proportional feeding.
 *
 * Evaluates how much resources cell has gained based on the cell's genotype
 * quality.
 * @param gain_from_the_genotype_ratio - the maximal permitted uptake rate, factor
 * dependent on the number of cells in the population used for scaling resources
 * available per
 * @param SharedBetweenFactor - factor dependent on the number of cells in the
 * population used for scaling resources available per one cell.
 *
 */
void Cell::evaluateInternalResourcesGainProportional(double gain_from_the_genotype_ratio,
        double SharedBetweenFactor) {
    double recentGainTemp;
    // gain for a cell in one iteration
    recentGainTemp = SharedBetweenFactor * gain_from_the_genotype_ratio \
        * CellsGrossGrowthRateAtTheMoment;
    recentGain = recentGainTemp;
    CellsResource += recentGainTemp;
    TotalGainSinceRepr += recentGainTemp;
}

/**
 * @brief Core method. Evaluates cell's internal resource cost.
 *
 * Evaluates how much resources cell has lost based on the cell's genotype
 * size and fixed parameters.
 * @param gene_replication_cost - cost of having one metabolic gene
 * @param core_genes - number of non-reducible genes
 * @param metabolic_costs - constant metabolic cost per cell
 */
void Cell::evaluateInternalResourcesCost(double gene_replication_cost,
        int core_genes, double metabolic_costs) {
    double GenomeSize = (double) CellsGenotype.GenotypeValues.size();
    double recentCostTemp;
    // sum of cost of living in one iteration
    recentCostTemp = -gene_replication_cost * ((double) core_genes \
        + GenomeSize + GenomeSize * GenomeSize) - metabolic_costs;
    if (recentCostTemp < -CellsResource) {
        recentLoss = -CellsResource;
        CellsResource = 0.0;
    }
    if (recentCostTemp >= -CellsResource) {
        recentLoss = recentCostTemp;
        CellsResource = CellsResource + recentCostTemp;
    }
}

/**
 * @brief Data collecting method. Sets the age of a new born cell.
 *
 * When called sets the age of the cell to ONE.
 */
void Cell::setCellAgeToOne() {
    CellsAge = 1.0;
}

/**
 * @brief Data collecting method. Increments the age of a cell.
 *
 * 'Ageing' function. When called adds one to cell's age.
 */
void Cell::ageingACell() {
    CellsAge += 1.0;
}

/**
 * @brief Data collecting method. Zero recent gain of a cell.
 *
 * Zeros the value of the cell's recent gain.
 */
void Cell::zeroRecentGainOfCell() {
    recentGain = 0.0;
}

/**
 * @brief Data collecting method. Unmarks a cell which ones reproduced.
 *
 * Changes the reproduction flag to 'false'.
 */
void Cell::setCellReprodFlagFalse() {
    DidReproduced = false;
    //    AvarageTotalGainTillRepr = 0.0;
    //    AvarageIntakeTillRepr = 0.0;
}

/**
 * @brief Data collecting method. Marks a cell which reproduced.
 *
 * Changes the reproduction flag to 'true'.
 */
void Cell::setCellReprodFlagTrue() {
    DidReproduced = true;
}

/**
 * @brief Data collecting method. Zeros the Cell's age passed since the
 * last reproduction.
 *
 * Zeros the Cell's age passed since the last reproduction.
 */
void Cell::setCellAARToZero() {
    AgeAtLastReprod = 0;
}

/**
 * @brief Data collecting method. Sets cell's resource at birth equal to the cell
 * resource variable.
 *
 * If this is a newborn cell, sets its resource at birth equal to the cell
 * resource variable.
 */
void Cell::setCellsResAfterReprod() {
    CellsResourceAfterReprod = CellsResource;
    TotalGainSinceRepr = 0.0;
}

/**
 * @brief Data collecting method. Zeros the average intake until reproduction.
 *
 * Zeros cell's average intake until reproduction.
 */
void Cell::zeroAvarageIntakeTillRepr() {
    AvarageIntakeTillRepr = 0.0;
}

/**
 * @brief Data collecting method. Zeros the average gain until reproduction.
 *
 * Zeros cell's average gain until reproduction.
 */
void Cell::zeroAvarageGainTillRepr() {
    AvarageTotalGainTillRepr = 0.0;
}

/**
 * @brief Core method / Data collecting method. Zeros all the parameter which
 * need to be zeroed for a new born cell
 *
 * Zeros all the parameter which need to be zeroed for a new born cell. These are:
 * TotalLoss, TotalGain, AvarageIntakeTillRepr, AgeSinceLastDivision,
 * AgeSinceLastDivision.
 */
void Cell::resetParamsOfNewCell() {
    CellsAge = 0.0;
    TotalLoss = 0.0;
    TotalGainSinceRepr = 0.0;
    AvarageIntakeTillRepr = 0.0;
    AvarageTotalGainTillRepr = 0.0;
    AgeAtLastReprod = 0;
    DidReproduced = false;
    NumberOfReprod = 0;
    CellsResourceAfterReprod = CellsResource;
}

/**
 * @brief Data collecting method. Sets the AgeAtLastReproduction variable to
 * the current age.
 *
 * Used to set the age of the cell for calculating resource circulation statistics
 */
void Cell::setAgeSinceLastReprod() {
    TimeSinceLastReprod = CellsAge - AgeAtLastReprod;
    AgeAtLastReprod = CellsAge;
}

/**
 * @brief Data collecting method. Sets the average net intake rate in the last period
 * preceding the current reproduction event.
 *
 * Sets the average net intake rate in the last period preceding the current
 * reproduction event.
 */
void Cell::setAvarageIntakeTillRepr() {
    AvarageIntakeTillRepr = (CellsResource - CellsResourceAfterReprod) / TimeSinceLastReprod;
}

/**
 * @brief Data collecting method. Sets the average total intake rate in the last period
 * preceding the current reproduction event.
 *
 * Sets the average real intake rate in the last period preceding the current
 * reproduction event.
 */
void Cell::setAvarageTotalGainTillRepr() {
    AvarageTotalGainTillRepr = TotalGainSinceRepr / TimeSinceLastReprod;
}

/**
 * @brief Data collecting method. Increments the number of reproductions of the cell.
 *
 * When called adds one to cell's number of reproductions.
 */
void Cell::incrementNumberOfReprod() {
    NumberOfReprod += 1;
}
