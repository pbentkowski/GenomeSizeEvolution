//      cell.h
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
#include <iostream>
#include <vector>

#ifndef CELL_H_
#define CELL_H_

/**
 * @class Cell
 * 
 * @brief Represents one cell. This is the agent of the model.
 *
 * It has a Genotype (with object of Gene class inside), stores resources, has
 * individual tag and tag of its mother. Has methods to grow, mutate, age. Also
 * some methods to allow to track the data record.
 */
class Cell {
public:
    // Core structures:
    // Cell's internal resource
    double CellsResource;
    // Cell's recent gain of resource
    double recentGain;
    // Cell's recent loss of resource
    double recentLoss;
    // Cell's genotype
    Genotype CellsGenotype;
    // number of iterations (time) since the cell's 'birth'
    double CellsAge;
    // the cell's growth rate at a current point of time
    double CellsGrossGrowthRateAtTheMoment;
    // did the cell mutated during last cycle, works bit like a flag
    bool TestIfMutated;
    
    // Data collecting structures:
    // tag which helps to calculate the speed of evolution
    unsigned long int OwnTag;
    // another tag which helps to calculate the speed of evolution
    unsigned long MothersTag;
    // yet another tag which helps to calculate the speed of evolution
    bool TestIfCounted;
    // tag helping to resolve reproduction questions
    bool DidReproduced;
    // number of times a cell reproduced
    int NumberOfReprod;
    // cell's age at the last division
    double AgeAtLastReprod;
    // time passed since the last division
    double TimeSinceLastReprod;
    // average intake of resource for a cell on its way to reproduce
    double AvarageIntakeTillRepr;
    // average total gain for a cell on its way to reproduce
    double AvarageTotalGainTillRepr;
    // Cell's resource at birth time
    double CellsResourceAfterReprod;
    // Resource flow data fields
    double TotalGainSinceRepr;
    double TotalLoss;

    Cell();
    // sets a new cell
    void setBrandNewCell(double Env, double minimal_cells_resource,
            double maximal_resourse_bonus_at_start, int min_genotype_initialized,
            int max_genotype_initialized, double not_super_gene_param,
            double gene_space_width, double environment_resolution,
            double one_over_total_env_space_surf);
    // decides (randomly) if cell mutates or not
    // if yes, calls mutation method from Genotype class
    void mutateCell(double Env, double duplication_probability,
            double deletion_probability, double mutate_threshold,
            int minimum_genotype_size, double not_super_gene_param,
            double gene_space_width, double environment_resolution,
            double one_over_total_env_space_surf);
    // calculates Gross Growth Rate for all the Env space
    void evaluateAdaptiveLandcape(double EnvArray[], int EnvSize);
    // calculates Gross Growth Rate just for one single Env value
    void evaluateGrowthRateFromGenotypeOnePointEnv(double Env);
    void setCellAgeToOne();
    void ageingACell();
    void zeroRecentGainOfCell();
    void setCellReprodFlagFalse();
    void setCellReprodFlagTrue();
    void setCellAARToZero();
    // calculates, basing on genotype quality and other costs, the gain of
    // internal resource when only few cells are randomly chosen to be
    // fed in a random sequence.
    void evaluateInternalResourcesGainSequential(double gain_from_the_genotype_ratio);
    // calculates, basing on genotype quality and other costs, the gain of
    // internal resource when all cells are fed in proportion to resource availability
    void evaluateInternalResourcesGainProportional(double gain_from_the_genotype_ratio,
            double InherentFraction);
    // calculates just costs
    void evaluateInternalResourcesCost(double gene_replication_cost, int core_genes,
            double metabolic_costs);
    void setCellsResAfterReprod();
    void zeroAvarageIntakeTillRepr();
    void zeroAvarageGainTillRepr();
    void resetParamsOfNewCell();
    void setAgeSinceLastReprod();
    void setAvarageIntakeTillRepr();
    void setAvarageTotalGainTillRepr();
    void incrementNumberOfReprod();
};

#endif
