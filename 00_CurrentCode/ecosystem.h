//      ecosystem.h
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
#include <iostream>
#include <vector>
#include <gsl/gsl_randist.h>

#ifndef ECOSYSTEM_H
#define ECOSYSTEM_H

/**
 * @struct GeneticVariability
 * 
 * @brief A struct which stores information for calculating the speed of evolution
 *
 * Stores two variables: a OwnTag and a flag showing was it already counted in an
 * appropriate function
 */
struct GeneticVariability {
    unsigned long int TheTag;
    bool IfCounted;
};

/**
 * @class Ecosystem
 * 
 * @brief This is 'The Ecosystem'. A place where things are happening.
 *
 * Stores a vector
 * (<a href="http://www.cplusplus.com/reference/stl/vector/">STL container</a>)
 * of cell objects. Has methods to set up a population of cells and later feed,
 * kill and reproduce them. Has a lot methods to harvest data.
 */
class Ecosystem {
public:
    // ----- Core fields -----
    std::vector<Cell> PopulationOfCells;
    double EnvFreeResources;

    // ----- Data collecting fields -----
    double meanGenomeSize;
    double STD_GenomeSize;
    unsigned long numberOfDeaths;
    unsigned long numberOfBirths;
    unsigned long startOfMutationCounter;

    // ----- Core methods -----
    Ecosystem();
    void setPopulation(double total_resource_pool, double Env,
            double minimal_cells_resource, double maximal_resourse_bonus_at_start,
            int min_genotype_initialized, int max_genotype_initialized,
            double not_super_gene_param, double gene_space_width,
            double environment_resolution, double one_over_total_env_space_surf);
    void reproduceOrKillCell(std::string mode_of_feeding, double Env,
            double duplication_probability, double deletion_probability,
            double mutate_threshold, int minimum_genotype_size,
            double not_super_gene_param, double gene_space_width,
            double environment_resolution, double one_over_total_env_space_surf,
            double random_death_rate, double minimal_cells_resource,
            double reproduction_resourse_size, double disparity_at_reproduction);
    void feedCells(std::string mode_of_feeding, double Env,
            double gene_replication_cost, int core_genes, double metabolic_costs,
            double gain_from_the_genotype_ratio);
    void updateGrowthRatesOfCellsInPopulation(double Env);
    void horizontalGeneTransfer(double hgt_donor_prob, double hgt_gene_trans);

    // ----- Data collecting methods -----
    void printNumberOfGenesInCells();
    double returnMeanNumberOfGenesInCells();
    double returnStdOfGenesInCells();
    double returnMeanAgeOfCells();
    double returnSTDAgeOfCells();
    double returnsSummOfCellResources();
    void histGenomeSize(int genome_size_histogram_size);
    void histCellsAge(int histSize);
    double returnPopulationShannonIndex();
    double returnGenesShannonIndex();
    void histGenotypeEnvelopeSurface(int envelope_of_genotype_hist_size,
            double envelope_of_genotype_hist_range);
    void histGenotypeEnvelopeMaximum(int maximum_of_genotype_hist_size,
            double maximum_of_genotype_hist_range);
    void histResourceAllocatedInCells(int resource_alloc_hist_size,
            double resource_alloc_hist_range, double resource_alloc_hist_start);
    void histResourceUptakenByCells(int resource_uptake_hist_size,
            double resource_uptake_hist_range);
    void histsAgeAndResourceAtReprod(int at_repr_intake_hist_size,
            double at_repr_intake_bin_range, int AAR_hist_size, double AAR_hist_range,
            int Reproduce_hist_size, double Reproduce_hist_range,
            int genome_size_histogram_size);
    double returnMeanEnvelopeToEnvRatio();
    double returnSTDEnvelopeToEnvRatio();
    void writeGeneralDataToFile(int k, double env);
    void cumulativeMutationNumbers(int current_time_stamp);
    void clearMutationCounters(int current_iteration);
    void writeGenotypesToFile(double ENV);
    void averageGenomeShape(double environment_resolution);
};

#endif /* ECOSYSTEM_H */
