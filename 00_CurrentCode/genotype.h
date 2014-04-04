//      genotype.h
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



#include "gene.h"
#include <iostream>
#include <vector>

#ifndef _GENOTYPE_H_
#define _GENOTYPE_H_

/**
 * @class Genotype
 *
 * @brief This is the genotype of one cell. Stores a vector
 * (<a href="http://www.cplusplus.com/reference/stl/vector/">STL container</a>)
 * of gene objects. Has methods to duplicate, delete and mutate them. Also calculates
 * growth rate for given environmental conditions.
 */
class Genotype {
public:
    // stores the genotype itself
    std::vector<Gene> GenotypeValues;
    // stores the shape of the genome's envelope (the U_i(x) function)
    std::vector<double> GenomeLandscape;
    // surface of Env space occupied by the genotype (envelope/total env ratio)
    double EnvSpaceFraction;
    // maximal value of genotype envelope
    double MaxGrowthRate;
    // total number of mutations from the beginning of all time which arose
    // in all of this genotype history
    unsigned int TotalPoinMutationNumber;
    unsigned int TotalDeletionNumber;
    unsigned int TotalDuplicationNumber;
    unsigned int TotalHGT_number;

    Genotype();
    void setGenotype(int min_genotype_initialized, int max_genotype_initialized,
            double not_super_gene_param, double gene_space_width,
            double environment_resolution);
    void deleteGenes(double deletion_probability, int minimum_genotype_size);
    void duplicateGene(double duplication_probability);
    void mutateGenes(double mutate_threshold, double not_super_gene_param,
            double gene_space_width);
    double calculateGrowthFromOneGene(int GenotypeValuesIndex, double Env);
    void evaluateAdaptiveLandcape(double environment_resolution,
            double one_over_total_env_space_surf);
    // data collecting methods
    void printGenotype();
    double getJustOneValueOfAGene(unsigned int GenotypeValuesIndex,
            unsigned int GeneValueIndex);
    void printNumberOfGenes();
    double returnNumberOfGenes();
    std::string genotypeToString();
};

#endif
