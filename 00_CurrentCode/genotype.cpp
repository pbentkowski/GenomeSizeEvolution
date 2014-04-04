//      genotype.cpp
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


#include "genotype.h"
#include "rngEngine.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <gsl/gsl_vector.h>

/**
 * @brief Core method. Constructor - sets a genotype, but with no values.
 * 
 * Just a constructor.
 * 	
 */
Genotype::Genotype() {
}

/**
 * @brief Core method. Sets up a brand new genotype. And assigns random number of
 * genes and values of the genes to it.
 *
 * Being given minimal and maximal permitted number of genes sets a genotype
 * and using methods from class Gene sets the genes values.
 *
 * @param min_genotype_initialized - minimal genotype size at init
 * @param max_genotype_initialized - maximal genotype at init
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - constant helping to calculate growth rate
 *
 */
void Genotype::setGenotype(int min_genotype_initialized, int max_genotype_initialized,
        double not_super_gene_param, double gene_space_width,
        double environment_resolution) {
    // calling and initialising random number generator
    // assigning random number of genes between 0 and 10
    rngEngine* p_rngEngine = rngEngine::getInstance();
    int Sizze;
    Gene tempGene;
    long int xSize = (long int) (2.0 / environment_resolution);
    do {
        Sizze = (int) ceil(min_genotype_initialized + (max_genotype_initialized \
            - min_genotype_initialized) * p_rngEngine->rngMTgetDauble());
    } while (Sizze == 0);
    // assigning values of genes
    for (int i = 0; i < Sizze; i++) {
        GenotypeValues.push_back(tempGene);
        GenotypeValues.back().setGeneValues(not_super_gene_param, gene_space_width);
    }
    GenomeLandscape.assign(xSize , 0.0); //add this to function's description 
    TotalPoinMutationNumber = 0;
    TotalDeletionNumber = 0;
    TotalDuplicationNumber = 0;
    TotalHGT_number = 0;
}

/**
 * @brief When called deletes some genes (at least one).
 *
 * Being given fixed deletion probability, randomly removes some genes,
 * but not less them one. This method is being called when The Cell 'decides'
 * it will delete a gene. When called it must delete at least one gene
 *
 * @param deletion_probability - gene deletion probability
 * @param minimum_genotype_size - minimal permitted size of the genotype
 *
 */
void Genotype::deleteGenes(double deletion_probability, int minimum_genotype_size) {
    // calling and initialising random number generator
    rngEngine* p_rngEngine = rngEngine::getInstance();
    //measuring the size of genotype
    int Size = (int) GenotypeValues.size();
    // marker of gene deletion
    bool TEST_IF_DELETED = false;
    // counter prevents loosing all the genes (and segmentation faults as a result)
    int counter = Size - minimum_genotype_size - 1;
    // calculating what's the probability of deletion of particular gene
    for (int i = Size - 1; i >= 0; i--) {
        if ((p_rngEngine->rngMTgetDauble() < deletion_probability) && (counter >= 0)) {
            GenotypeValues.erase(GenotypeValues.begin() + i);
            TEST_IF_DELETED = true;
            counter -= 1;
            TotalDeletionNumber += 1;
        }
    }
    // if TEST_IF_DELETED is false it means that there were no deletions,
    // so the first gene is being deleted by force
    if ((TEST_IF_DELETED == false) && (counter >= 0)) {
        GenotypeValues.erase(GenotypeValues.begin());
        TotalDeletionNumber += 1;
    }
    // if TEST_IF_DELETED is true it means that there were deletions already,
    // so a dice can be thrown to decide if the first gene will be deleted
    if ((TEST_IF_DELETED == true) && (p_rngEngine->rngMTgetDauble() \
            < deletion_probability) && (counter >= 0)) {
        GenotypeValues.erase(GenotypeValues.begin());
        TotalDeletionNumber += 1;
    }
    TEST_IF_DELETED = false;
}

/**
 * @brief When called duplicates some genes (at least one duplicate ).
 *
 * Being given fixed duplication probability from file defines.h randomly
 * duplicates some genes, but not less them one. This method is being called
 * when The Cell 'decides' it will duplicate a gene. When called it must duplicate
 * at least one gene
 * @param duplication_probability - gene duplication probilty
 *
 */
void Genotype::duplicateGene(double duplication_probability) {
    // calling and initialising random number generator
    rngEngine* p_rngEngine = rngEngine::getInstance();
    // measuring the size of genotype
    int Size = (int) GenotypeValues.size();
    // marker of gene duplication
    bool TEST_IF_DUPLICATED = false;
    // calculating what's the probability of duplication of particular gene
    for (int i = 0; i < Size - 1; i++) {
        if (p_rngEngine->rngMTgetDauble() < duplication_probability) {
            GenotypeValues.push_back(GenotypeValues[i]);
            TEST_IF_DUPLICATED = true;
            TotalDuplicationNumber += 1;
        }
    }
    // if TEST_IF_DUPLICATED is false it means that there were no duplications,
    // so the last gene is being duplicated by force
    if (TEST_IF_DUPLICATED == false) {
        GenotypeValues.push_back(GenotypeValues.back());
        TotalDuplicationNumber += 1;
    }
    // if TEST_IF_DUPLICATED is true it means that there were duplications
    // already, so a dice can be thrown to decide if last gene will be duplicated
    if ((TEST_IF_DUPLICATED == true) && (p_rngEngine->rngMTgetDauble() \
            < duplication_probability)) {
        GenotypeValues.push_back(GenotypeValues.back());
        TotalDuplicationNumber += 1;
    }
    TEST_IF_DUPLICATED = false;
}

/**
 * @brief Core method. When called mutates some genes (at least one).
 *
 * Being given fixed mutation probability, randomly mutates
 * some genes using mutation method from class Gene, but not less them one in
 * each operation. This method is being called when The Cell 'decides' it
 * will mutate a gene. When called it must mutate at least one gene.
 * 
 * @param mutate_threshold - gene modification probalitily
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - he width of the env space genes can take
 *
 */
void Genotype::mutateGenes(double mutate_threshold, double not_super_gene_param,
        double gene_space_width) {
    // calling and initializing random number generator
    rngEngine* p_rngEngine = rngEngine::getInstance();
    // measuring the size of genotype
    int Size = (int) GenotypeValues.size();
    // marker of gene mutation
    bool TEST_IF_MUTATED = false;
    // calculating what's the probability of mutating of particular gene
    for (int i = 0; i < Size - 1; i++) {
        if (p_rngEngine->rngMTgetDauble() < mutate_threshold) {
            GenotypeValues[i].setGeneValues(not_super_gene_param, gene_space_width);
            TEST_IF_MUTATED = true;
            TotalPoinMutationNumber += 1;
        }
    }
    // if TEST_IF_MUTATED is false it means that there were no mutations,
    // so the last gene is being mutated by force
    if (TEST_IF_MUTATED == false) {
        GenotypeValues.back().setGeneValues(not_super_gene_param, gene_space_width);
        TotalPoinMutationNumber += 1;
    }
    // if TEST_IF_MUTATED is true it means that there were mutations already,
    // so a dice can be thrown to decide if last gene will be mutated
    if ((TEST_IF_MUTATED == true) && (p_rngEngine->rngMTgetDauble() < mutate_threshold)) {
        GenotypeValues.back().setGeneValues(not_super_gene_param, gene_space_width);
        TotalPoinMutationNumber += 1;
    }
    TEST_IF_MUTATED = false;
}

/**
 * @brief Core method. Calculates growth rate for given genes in given conditions.
 *
 * Calculate growth rate from one gene, which index in Genotype vector is
 * given as an argument for environmental conditions given as a second argument.
 * Uses method from class Gene.
 * @param GenotypeValuesIndex - position of a gene in genotype, 
 * @param Env - current value of the environmental conditions.
 * @return GenotypeValues or -1 (error message)
 *
 */
double Genotype::calculateGrowthFromOneGene(int GenotypeValuesIndex, double Env) {
    if (((int) GenotypeValues.size()) <= GenotypeValuesIndex) {
        std::cout << "GenotypeValuesIndex provided in calculateGrowthFromOneGene " << \
            "(int, double) is outside GenotypeValues vector range" << std::endl;
        return -1;
    } else {
        return GenotypeValues[GenotypeValuesIndex].calculateGrowthFromGene(Env);
    }
}

/**
 * @brief Core method. Returns one given value in gene array of a gene in a
 * given position.
 * 
 * @param GenotypeValuesIndex - position of a gene in genotype
 * @param GeneValueIndex - position of one specified value out of three a gene has
 * @return requested value (double) or -1 (error message)
 */
double Genotype::getJustOneValueOfAGene(unsigned int GenotypeValuesIndex,
        unsigned int GeneValueIndex) {
    if (((unsigned int) GenotypeValues.size()) <= GenotypeValuesIndex) {
        std::cout << "GenotypeValuesIndex provided in getJustOneValueOfAGene" << \
            "(int, int) is outside GenotypeValues vector range" << std::endl;
        return -1;
    } else {
        return GenotypeValues[GenotypeValuesIndex].getOneValue(GeneValueIndex);
    }
}


/**
 * @brief Test method. Print the genotype to the screen.
 *
 * Test function. Not used in 'regular' run.
 */
void Genotype::printGenotype() {
    std::cout << "% number of genes in genotype: " << \
        (int) GenotypeValues.size() << std::endl;
    for (int i = 0; i < (int) GenotypeValues.size(); i++) {
        std::cout << GenotypeValues[i].geneToString() << std::endl;
    }
}

/**
 * @brief Data collecting method. Puts a number of genes in the genotype.
 *
 */
void Genotype::printNumberOfGenes() {
    std::cout/* << "% number of genes in genotype: "*/ << GenotypeValues.size() \
        << " " /* std::endl*/;
}

/**
 * @brief Data collecting method. Returns the number of genes in the genotype.
 *
 */
double Genotype::returnNumberOfGenes() {
    return (double) GenotypeValues.size();
}

/** \brief Data collecting method. Counts the coverage of environmental
 * space surface covered by the surface of the envelope of genotype.
 *
 * Calculates how much of the total environment space is being covered by
 * surface under genotype's envelope. Uses <a href=
 * "http://www.gnu.org/software/gsl/manual/html_node/Vectors.html">
 * GNU Scientific Library vectors</a> to compute those values.
 * 
 * @param environment_resolution - resolution of the environment space
 * @param one_over_total_env_space_surf - 1/total env surface (const supporting computation).
 */
void Genotype::evaluateAdaptiveLandcape(double environment_resolution,
        double one_over_total_env_space_surf) {
    // measuring genome size and env resolution
    int GenomeSize = (int) GenotypeValues.size();
    long int GLsize = (long int) (2.0 / environment_resolution);
    // env value
    double env = -1.0;
    // surface
    double s = 0.0;
    MaxGrowthRate = 0.0;
    double maxGeneValue;
    // Temporary vector collecting Gross Growth Rates provided by all the
    // genes for a given Environment conditions
    gsl_vector * TempVector = gsl_vector_alloc(GenomeSize);
    for (int e = 0; e < GLsize; e++) {       // while (env <= 1.0)
        for (int j = 0; j < GenomeSize; j++) {
            // calculating growth rate for a given gene in a given Environment value
            gsl_vector_set(TempVector, j, calculateGrowthFromOneGene(j, env));
        }
        maxGeneValue = gsl_vector_max(TempVector);
        GenomeLandscape[e] = maxGeneValue;
        if (maxGeneValue > MaxGrowthRate) {
            // finding the global maximum of the genotype's envelope
            MaxGrowthRate = maxGeneValue;
        }
        // calculating the surface of one bin using highest gene value and
        // incrementing total surface
        s = s + maxGeneValue * environment_resolution;
        env = env + environment_resolution;
    }
    EnvSpaceFraction = s * one_over_total_env_space_surf;
    gsl_vector_free(TempVector);
}

/**
 * @brief Data collecting support method. Transforms the list of genotype's genes
 * into a <a href="http://www.cplusplus.com/reference/string/string/">STD string</a>.
 * 
 * @return all three values of all the genes in the genotype arranged in a 
 * <a href="http://www.cplusplus.com/reference/string/string/">STD string</a>.
 */
std::string Genotype::genotypeToString(){
    std::ostringstream genes;
    int Size = (int) GenotypeValues.size();
    for (int i = 0; i < Size; i++) {
        genes << GenotypeValues[i].geneToString() << " ";
    }
    return genes.str();
}