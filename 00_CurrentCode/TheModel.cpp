//      TheModel.cpp
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

// Compilation:
// g++ -O2 -o Model_with_HGT -L/usr/local/lib TheModel.cpp genotype.cpp gene.cpp
// cell.cpp ecosystem.cpp rngEngine.cpp tagging_system.cpp -lgsl -lgslcblas -lm


#include "ecosystem.h"
#include "tagging_system.h"
#include <math.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "rngEngine.h"

/**
 * @brief  Function generating new values of the abiotic environment.
 * 
 * @param turbulence - how big is the maximal possible change
 * @param twice_turbulence - twice the value above
 * @param RNG - a number between [0, 1] randomly generated elsewhere 
 * @param env - the old value of the environmental conditions
 * @return env - a new value of the environmental conditions
 */
double envGenerate(double turbulence, double twice_turbulence, double RNG,
        double env){
    double delta_env;
    double tmp_env;
    delta_env = twice_turbulence * RNG - turbulence;
    tmp_env = env + delta_env;
    if (tmp_env < -1.0) {
        env = -2.0 - (env + delta_env);
    }
    else if (tmp_env > 1.0) {
        env = 2.0 - (env + delta_env);
    } else {
        env = tmp_env;
    }
    return env;
}

/**
 * @brief  The main function.
 * 
 * @param argc - number of parameters at input
 * @param argv - the good old char *argv[] (43 parameters)
 * @return  0 if succesfully run
 */
int main(int argc, char *argv[]) {

    time_t start, end;
    double dif;
    time(&start);

    if (argc < 45) {
        std::cout << "Not enough arguments. It has to be " << \
            "precisely 45 of them and only " << argc << " were provided." << std::endl;
        return 0;
    }
    if (argc > 45) {
        std::cout << "Too many arguments. It has to be " <<
                "precisely 45 of them but " << argc << " were provided." << std::endl;
        return 0;
    }

    // loading general parameters
    // a measure of stability of the env:
    double turbulence = atof(argv[1]);
    double twice_turbulence = turbulence * 2.0;
    // Environment has values [-1:1], this is its resolution:
    double environment_resolution = atof(argv[2]);
    // it's 1/(1*(1-(-1))):
    double one_over_total_env_space_surf = atof(argv[3]);
    // Random number generator seed :
    int rng_seed = atoi(argv[4]);
    // number of iterations of main loop:
    int number_of_iterations = atoi(argv[5]);
    // loading gene class parameters
    // parameter which scales the width of Gaussian curve to its hight
    // (prevents super-genes to appear):
    double not_super_gene_param = atof(argv[6]);
    // probability of mutation in single gene (change of the gene's
    // Gaussian curve shape):
    double mutate_threshold = atof(argv[7]);
    // how much wither is the gene space comparing to environmental space,
    // generally place 1:
    double gene_space_width = atof(argv[8]);
    // loading genotype class parameters
    // deletion_probability must be between [0:1]
    double deletion_probability = atof(argv[9]);
    // duplication_probability must be between [0:1]
    double duplication_probability = atof(argv[10]);
    // maximal number of genes at initialisation
    int max_genotype_initialized = atoi(argv[11]);
    // minimal number of genes at initialisation
    int min_genotype_initialized = atoi(argv[12]);

    //  minimum_genotype_size is the minimal number (integer) of genes
    // a genotype can have (basically it prevents segmentation faults) !!!
    // IT CANNOT BE SMALLER THAN 1 (one) !!! Smaller value may result in
    // segmentation faults:
    int minimum_genotype_size = atoi(argv[13]);
    // loading cell class parameters
    // number of core genes:
    int core_genes = atoi(argv[14]);
    // cost of maintaining and expressing one gene:
    double gene_replication_cost = atof(argv[15]);
    // how much cell can gain from having an maximal 'fitness' at particular
    // Env conditions:
    double gain_from_the_genotype_ratio = atof(argv[16]);
    // constant cost of living:
    double metabolic_costs = atof(argv[17]);
    // below this value cell dies:
    double minimal_cells_resource = atof(argv[18]);
    // at this size cell duplicates:
    double reproduction_resourse_size = atof(argv[19]);
    // maximal value of randomly assigned resources overhead when a cell
    // is initialised at the beginning of the model:
    double maximal_resourse_bonus_at_start = atof(argv[20]);
    // loading ecosystem class parameters
    // dying because of viruses, predation, bad weather, cosmic radiation etc.:
    double random_death_rate = atof(argv[21]);
    // how big is the resource poll (limiting factor for population size):
    double total_resource_pool = atof(argv[22]);
    // loading data harvesting parameters
    // size of the histogram of gnome sizes, has to be as big as predicted
    // number of genes but not
    // smaller then MAX_GENOTYPE_INITIALIZED:
    int genome_size_histogram_size = atoi(argv[23]);
    // size of the histogram of envelope of genotype to environmental space ratio:
    int envelope_of_genotype_hist_size = atoi(argv[24]);
    // range of bins of the histogram of envelope of genotype to
    // environmental space ratio:
    double envelope_of_genotype_hist_range = atof(argv[25]);
    // size of the histogram of maximums of genotypes' envelopes:
    int maximum_of_genotype_hist_size = atoi(argv[26]);
    // range of bins of the histogram of maximums of genotypes' envelopes:
    double maximum_of_genotype_hist_range = atof(argv[27]);
    // size of the histogram of ages of cells:
    int AgeHistSize = atoi(argv[28]);
    // size of the histogram of resource allocated to cells:
    int resource_alloc_hist_size = atoi(argv[29]);
    // width of bins of the histogram of resource allocated to cells:
    double resource_alloc_hist_range = atof(argv[30]);
    // beginning value of the histogram of resource allocated to cells:
    double resource_alloc_hist_start = atof(argv[31]);
    // size of the histogram of resource taken by cells:
    int resource_uptake_hist_size = atoi(argv[32]);
    // width of bins of the histogram of resource taken by cells:
    double resource_uptake_hist_range = atof(argv[33]);
    // sets the mode of feeding:
    std::string mode_of_feeding(argv[34]);
    // how uneven is resource distribution when cell is dividing:
    double disparity_at_reproduction = atof(argv[35]);
    // size of the histogram of the age at the first reproduction:
    int AAR_hist_size = atoi(argv[36]);
    // width of bins of the histogram of the age at the first reproduction:
    double AAR_hist_range = atof(argv[37]);
    // size of the histogram of the number of reproduction events:
    int Reproduce_hist_size = atoi(argv[38]);
    // width of bins of the histogram of the number of reproduction events:
    double Reproduce_hist_range = atof(argv[39]);
    // size of the histogram of the average intake per unit of time until
    // reproduction:
    int at_repr_intake_hist_size = atoi(argv[40]);
    // width of bins of the histogram  of the average intake per unit of
    // time until reproduction:
    double at_repr_intake_bin_range = atof(argv[41]);
    // Horizontal gene transfer parameters:
    double hgt_donor_prob = atof(argv[42]);
    double hgt_gene_trans = atof(argv[43]);
    // Interval of data sampling for the output files
    int sampling_interv = atoi(argv[44]);
    
    // testing if entered parameters make sense
    bool param_does_not_make_sense_flag = false;
    if (turbulence > 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Value of turbulence level (argument #1) larger then" << \
            " 1.0 makes no sense. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (environment_resolution >= 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Value of environment's resolution (argument #2) " << \
            "cannot be larger or equal 1.0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (one_over_total_env_space_surf != 0.5) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Value of one_over_total_env_space_surf " << \
            "(argument #3) has to be 0.5. Sorry. Check entered value." \
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (number_of_iterations <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "NUmber of iterations (argument #5) has to be an " << \
            "integer larger then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (not_super_gene_param <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "'Non super gene' parameter (argument #6) has to " << \
            "be larger then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (mutate_threshold < 0.0 || mutate_threshold > 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Probability of single gene mutation (argument #7) " << \
            "has to be between 0 and 1. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (gene_space_width <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Gene space width (argument #8) has to be larger " << \
            "then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (deletion_probability < 0.0 || deletion_probability > 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Probability of deleting a single gene " << \
            "(argument #9) has to be between 0 and 1. Check entered value." \
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (duplication_probability < 0.0 || duplication_probability > 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Probability of duplication of a single gene " << \
            "(argument #10) has to be between 0 and 1. Check entered value."\
        << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (max_genotype_initialized <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Maximal number of gene at start(argument #11) " << \
            "has to be larger then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (min_genotype_initialized <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Minimal number of gene at start(argument #11) " << \
            "has to be larger then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (min_genotype_initialized > max_genotype_initialized) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Minimal number of genes at initialisation " << \
            "(argument #12) cannot be larger than maximal number of " << \
            "genes at initialisation (parameter #11). Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (minimum_genotype_size <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Minimal number of genes at all times " << \
            "(argument #13) has to be larger then 0. Check entered value."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (core_genes < 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Number of core genes (argument #14) cannot be " << \
            "smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (gene_replication_cost <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The cost of maintaining a gene (argument #15) " << \
            "cannot be smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (gain_from_the_genotype_ratio <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The gain from genotype ratio (argument #16) " << \
            "cannot be smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (metabolic_costs <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The metabolic costs (argument #17) cannot be " << \
            "smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (minimal_cells_resource <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Minimal cell resource to live (argument #18) " << \
            "cannot be smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (reproduction_resourse_size <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Minimal cell resource to reproduce (argument #19) " << \
            "cannot be smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (reproduction_resourse_size < minimal_cells_resource) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Minimal cell resource to reproduce (argument #19) " << \
            "cannot be smaller then minimal cell resource to live " << \
            "(argument #18). Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (maximal_resourse_bonus_at_start < 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Maximal resource bonus at start (argument #20) " << \
            "cannot be smaller then 0. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (random_death_rate < 0.0 || random_death_rate > 1) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Probability of random death (argument #21) has " << \
            "to be between 0 and 1. Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (total_resource_pool < 0.0 || total_resource_pool < reproduction_resourse_size \
        || total_resource_pool < minimal_cells_resource) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "Total amount of resource in the system " << \
            "(argument #22) has to be larger then 0, larger then minimal " << \
            "cell resource to reproduce and larger then minimal cell " << \
            "resource to live. Well, it has to be quite big in general. " << \
            "Check entered value." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (genome_size_histogram_size < max_genotype_initialized) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the histogram of gnome  " << \
            "sizes (argument #23) cannot be smaller than maximal number  " << \
            "of genes at initialisation (argument #11). Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (resource_alloc_hist_size <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the histogram of resource " << \
            "allocated in the cells (argument #29) cannot be smaller " << \
            "than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (resource_alloc_hist_range <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the bin of histogram of " << \
            "resource allocated in the cells (argument #30) cannot be " << \
            "smaller than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (resource_alloc_hist_start <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of initial value of histogram of " << \
            "resource allocated in the cells (argument #30) cannot be " << \
            "smaller than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if ((mode_of_feeding.compare("p") != 0) && (mode_of_feeding.compare("s") != 0)) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The feeding mode (argument #34) has to be a string " << \
            "\'p\' (for proportional feeding) or \'s\' (for sequential " << \
            "feeding). Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (disparity_at_reproduction > 0.5 || disparity_at_reproduction < 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The disparity of resource at reproduction " << \
            "(argument #35) has to be larger then 0 and smaller " << \
            "then 0.5. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (AAR_hist_size <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the histogram of the age " << \
            "at the first reproduction (argument #36) cannot be smaller " << \
            "than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (AAR_hist_range <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the bin of histogram of " << \
            "the age at the first reproduction of the cells (argument #37) " << \
            "cannot be smaller than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (Reproduce_hist_size <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the histogram of the number " << \
            "of reproductions of the cells (argument #38) cannot be " << \
            "smaller than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (Reproduce_hist_range <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the bin of histogram of " << \
            "of the number of reproduction events the cells (argument #39) " << \
            "cannot be smaller than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (at_repr_intake_hist_size <= 0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of size of the histogram of the average " << \
            "intake per unit of time until reproduction (argument #40) " << \
            "cannot be smaller than 0. Check entered values." << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (at_repr_intake_bin_range <= 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of the bin width of the histogram of " << \
            "the average intake per unit of time until reproduction " << \
            "(argument #41) cannot be smaller than 0. Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if ( hgt_donor_prob < 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of the probability of being a donor in HGT" << \
            "(argument #42) cannot be smaller than 0. Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if ( hgt_donor_prob > 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of the probability of being a donor in HGT" << \
            "(argument #42) cannot be larger then 1. Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if ( hgt_gene_trans < 0.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of the probability of transfering a gene in " << \
            "HGT (argument #43) cannot be smaller than 0. Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if ( hgt_gene_trans > 1.0) {
        std::cout << "Error in defining parameters:" << std::endl;
        std::cout << "The value of the probability of transfering a gene in " << \
            "HGT (argument #43) cannot be larger then 1. Check entered values."\
            << std::endl;
        param_does_not_make_sense_flag = true;
    }
    if (param_does_not_make_sense_flag == true) {
        std::cout << "There were errors in defining parameters. Run aborted."\
            << std::endl;
        std::cout << "Use parameter set-up program to enter parameters."\
            << std::endl;
        return 0;
    }
    if (param_does_not_make_sense_flag == false) {
        std::cout << "Everything seems fine... Let's go!" << std::endl;
    }

    // initialising singleton classes (RNG and the tag generator)
    rngEngine* p_rngEngine = rngEngine::getInstance();
    if (rng_seed == -1) {
        p_rngEngine->rngMTset(time(0));
    } else {
        p_rngEngine->rngMTset(rng_seed);
    }
    Tagging_system* pTagging_system = Tagging_system::getInstance();
    pTagging_system->setValue(0);
    // initialising environmental conditions value
    double env = 0.0;
//    double tmp_env = 0.0;
//    double delta_env = 0.0;
    int half_numb_of_iterations = number_of_iterations / 2 ;

    //-----Writing parameters to file storing parameters -----------------------
    std::ofstream ModelParams;
    ModelParams.open("ModelParams.dat");
    ModelParams << "# This file stores informations about model's run parameters "\
        << std::endl;
    ModelParams << "# This model run aims in checking is how turbulence impact " << \
        "Shannon index.\n#enironment was given with a function: env(t+1) " << \
        "= env(t) + (turbulence*rand(0,1) - turbulence/2)" << std::endl;
    ModelParams << std::endl;
    ModelParams << "# general: " << std::endl;
    ModelParams << "    turbulence = " << turbulence << std::endl;
    ModelParams << "    environment_resolution = " << environment_resolution << std::endl;
    ModelParams << "    rng_seed = " << rng_seed << std::endl;
    ModelParams << "    number of iterations = " << number_of_iterations << std::endl;
    ModelParams << "# for gene class: " << std::endl;
    ModelParams << "    not_super_gene_param = " << not_super_gene_param << std::endl;
    ModelParams << "    mutate_threshold = " << mutate_threshold << std::endl;
    ModelParams << "    gene_space_width = " << gene_space_width << std::endl;
    ModelParams << "# for genotype class: " << std::endl;
    ModelParams << "    deletion_probability = " << deletion_probability << std::endl;
    ModelParams << "    duplication_probability = " << duplication_probability << std::endl;
    ModelParams << "    max_genotype_initialised = " << max_genotype_initialized << std::endl;
    ModelParams << "    min_genotype_initialised = " << min_genotype_initialized << std::endl;
    ModelParams << "    minimum_genotype_size = " << minimum_genotype_size << std::endl;
    ModelParams << "# for cell class: " << std::endl;
    ModelParams << "    core_genes = " << core_genes << std::endl;
    ModelParams << "    gene_replication_cost = " << gene_replication_cost << std::endl;
    ModelParams << "    gain_from_the_genotype_ratio = " << gain_from_the_genotype_ratio << std::endl;
    ModelParams << "    metabolic_costs = " << metabolic_costs << std::endl;
    ModelParams << "    minimal_cells_resource = " << minimal_cells_resource << std::endl;
    ModelParams << "    reproduction_resourse_size = " << reproduction_resourse_size << std::endl;
    ModelParams << "    maximal_resourse_bonus_at_start = " << maximal_resourse_bonus_at_start << std::endl;
    ModelParams << "# for ecosystem class: " << std::endl;
    ModelParams << "    random_death_rate = " << random_death_rate << std::endl;
    ModelParams << "    total_resource_pool = " << total_resource_pool << std::endl;
    ModelParams << "# data collecting const : " << std::endl;
    ModelParams << "    envelope_of_genotype_hist_size = " << envelope_of_genotype_hist_size << std::endl;
    ModelParams << "    envelope_of_genotype_hist_range = " << envelope_of_genotype_hist_range << std::endl;
    ModelParams << "    maximum_of_genotype_hist_size = " << maximum_of_genotype_hist_size << std::endl;
    ModelParams << "    maximum_of_genotype_hist_range = " << maximum_of_genotype_hist_range << std::endl;
    ModelParams << "    AgeHistSize = " << AgeHistSize << std::endl;
    ModelParams << "    resource_alloc_hist_size = " << resource_alloc_hist_size << std::endl;
    ModelParams << "    resource_alloc_hist_range = " << resource_alloc_hist_range << std::endl;
    ModelParams << "    resource_alloc_hist_start = " << resource_alloc_hist_start << std::endl;
    ModelParams << "    resource_uptake_hist_size = " << resource_uptake_hist_size << std::endl;
    ModelParams << "    resource_uptake_hist_range = " << resource_uptake_hist_range << std::endl;
    ModelParams << "# Mode of feeding (p - proportional, s - sequential)" << std::endl;
    ModelParams << "    mode_of_feeding = " << mode_of_feeding << std::endl;
    ModelParams << "    disparity_at_reproduction = " << disparity_at_reproduction << std::endl;
    ModelParams << "    AFR_hist_size = " << AAR_hist_size << std::endl;
    ModelParams << "    AFR_hist_range = " << AAR_hist_range << std::endl;
    ModelParams << "    Reproduce_hist_size = " << Reproduce_hist_size << std::endl;
    ModelParams << "    Reproduce_hist_range = " << Reproduce_hist_range << std::endl;
    ModelParams << "    at_repr_intake_hist_size = " << at_repr_intake_hist_size << std::endl;
    ModelParams << "    at_repr_intake_bin_range = " << at_repr_intake_bin_range << std::endl;
    ModelParams << "    genome_size_histogram_size = " << genome_size_histogram_size << std::endl;
    ModelParams << "# horizontal gene transfer params: " << std::endl;
    ModelParams << "    hgt_donor_prob = " << hgt_donor_prob << std::endl;
    ModelParams << "    hgt_gene_trans = " << hgt_gene_trans << std::endl;
    ModelParams.close();
    //------------------------------------------------------------

    // ------ Creating files with data ---------------------------
    std::ofstream GetotypeHistFile;
    GetotypeHistFile.open("GenomeSizeData.dat");
    GetotypeHistFile.close();

    std::ofstream GeneralDataFile;
    GeneralDataFile.open("GeneralData.dat");
    GeneralDataFile << "#time env_conditions mean_num_genes " << \
        "STD_of_num_genes num_of_cells res_in_cells res_in_Env " << \
        "mean_age STD_age Shannon_index born killed Envelope_mean " << \
        "Envelope_STD" << std::endl;
    GeneralDataFile.close();

    std::ofstream CellsAgeFile;
    CellsAgeFile.open("CellsAgeData.dat");
    CellsAgeFile << "# age of the cell - histograms" << std::endl;
    CellsAgeFile.close();

    std::ofstream FrameSizeFile;
    FrameSizeFile.open("FrameSizeData.dat");
    FrameSizeFile << "# size of the genotype's frame - histograms" << std::endl;
    FrameSizeFile.close();

    std::ofstream FrameMaxFile;
    FrameMaxFile.open("FrameMaxData.dat");
    FrameMaxFile << "# maximums of the genotype's frame - histograms"\
        << std::endl;
    FrameMaxFile.close();

    std::ofstream ResourceInCellsFile;
    ResourceInCellsFile.open("RecourceInCells.dat");
    ResourceInCellsFile << "# resources allocated into cells in one " << \
        "time step - histograms" << std::endl;
    ResourceInCellsFile.close();

    std::ofstream ResourceUptakeCellsFile;
    ResourceUptakeCellsFile.open("RecourceUptakenCells.dat");
    ResourceUptakeCellsFile << "# resources taken by cells in one time " << \
        "step - histograms" << std::endl;
    ResourceUptakeCellsFile.close();

    std::ofstream AgeAtFirstReproductionFile;
    AgeAtFirstReproductionFile.open("AgeTillReproduction.dat");
    AgeAtFirstReproductionFile << "# Time passed since the last time cell " << \
        "reproduced (or was born) - histograms" << std::endl;
    AgeAtFirstReproductionFile.close();

    std::ofstream NumberOfReproductionsFile;
    NumberOfReproductionsFile.open("NumberOfReproductions.dat");
    NumberOfReproductionsFile << "# How many times a cell reproduced " << \
        "- histograms" << std::endl;
    NumberOfReproductionsFile.close();

    std::ofstream AverageIntakeTillReprodFile;
    AverageIntakeTillReprodFile.open("AvarIntakeAtRepr.dat");
    AverageIntakeTillReprodFile << "# Cells' average intake before " << \
        "reproduction - histograms" << std::endl;
    AverageIntakeTillReprodFile.close();

    std::ofstream AverageTotGainTillReprodFile;
    AverageTotGainTillReprodFile.open("AvarGainAtRepr.dat");
    AverageTotGainTillReprodFile << "# Cells' average total gain before " << \
        "reproduction - histograms" << std::endl;
    AverageTotGainTillReprodFile.close();

    std::ofstream GeneNumberOfReprodFile;
    GeneNumberOfReprodFile.open("GeneNumberOfRepr.dat");
    GeneNumberOfReprodFile << "# Cells' gene number distribution. Only cells " << \
        "which reproduced - histograms" << std::endl;
    GeneNumberOfReprodFile.close();

    std::ofstream MutationRecords;
    MutationRecords.open("MutationRecords.dat");
    MutationRecords << "# cumulative number of mutations" << std::endl;
    MutationRecords << "# OwnTag PoinMutations Duplications Deletions " << \
        "HGT_accepted" << std::endl;
    MutationRecords.close();
    
    std::ofstream AllTheGenomes;
    // ooop, watch that line below when changing the file name!
    AllTheGenomes.open("GenomesOfPopulation.dat");
    AllTheGenomes << "# these are all the genomes in the ecosystem" << std::endl;
    AllTheGenomes << "#number_of_cells_in_clone mothers_tag own_tag gene1[0] "
            "gene1[1] gene1[2] ... geneN[0] geneN[1] geneN[2]" << std::endl;
    AllTheGenomes.close();
    
    std::ofstream histGenotypeFile;
    histGenotypeFile.open("AvaregeGenotype.dat");
    histGenotypeFile << "# Histogram of U(x) values (intake values) in all the"
            " cells in the population" << std::endl;
    histGenotypeFile.close();
    
    std::ofstream GenotEnvelMean;
    GenotEnvelMean.open("GenotypeEnvelMean.dat");
    GenotEnvelMean << "# Mean of U(x) function (the genome's envelope)" << std::endl;
    GenotEnvelMean.close();
    
    std::ofstream GenotEnvelSTD;
    GenotEnvelSTD.open("GenotypeEnvelSTD.dat");
    GenotEnvelSTD << "# STD of U(x) function (the genome's envelope)" << std::endl;
    GenotEnvelSTD.close();
    // -----------------------------------------------------

    // ----------- The main stuff -------------
    Ecosystem A_Ecosystem;
    A_Ecosystem.setPopulation(total_resource_pool, env, minimal_cells_resource,
            maximal_resourse_bonus_at_start, min_genotype_initialized,
            max_genotype_initialized, not_super_gene_param, gene_space_width,
            environment_resolution, one_over_total_env_space_surf);
    A_Ecosystem.writeGeneralDataToFile(0, env);  // better not switch off this one
    A_Ecosystem.histGenomeSize(genome_size_histogram_size);
    A_Ecosystem.histCellsAge(AgeHistSize);
    A_Ecosystem.histGenotypeEnvelopeSurface(envelope_of_genotype_hist_size,
            envelope_of_genotype_hist_range);
    A_Ecosystem.histGenotypeEnvelopeMaximum(maximum_of_genotype_hist_size,
            maximum_of_genotype_hist_range);
    A_Ecosystem.histResourceAllocatedInCells(resource_alloc_hist_size,
            resource_alloc_hist_range, resource_alloc_hist_start);
    A_Ecosystem.histResourceUptakenByCells(resource_uptake_hist_size,
            resource_uptake_hist_range);
    A_Ecosystem.histsAgeAndResourceAtReprod(at_repr_intake_hist_size,
            at_repr_intake_bin_range, AAR_hist_size, AAR_hist_range,
            Reproduce_hist_size, Reproduce_hist_range, genome_size_histogram_size);
    A_Ecosystem.averageGenomeShape(environment_resolution);

    for (int k = 1; k <= half_numb_of_iterations; k++) {
        // Tossing the next value of the environmental conditions
        env = envGenerate(turbulence, twice_turbulence, p_rngEngine->rngMTgetDauble(), env);
        // all the stuff related to living
        A_Ecosystem.feedCells(mode_of_feeding, env, gene_replication_cost, core_genes,
                metabolic_costs, gain_from_the_genotype_ratio);
        // data extraction from the simulation
        if (k % sampling_interv == 0) {
            A_Ecosystem.writeGeneralDataToFile(k, env);  // better not switch off this one
            A_Ecosystem.histGenomeSize(genome_size_histogram_size);
            A_Ecosystem.histCellsAge(AgeHistSize);
            A_Ecosystem.histGenotypeEnvelopeSurface(envelope_of_genotype_hist_size,
                    envelope_of_genotype_hist_range);
            A_Ecosystem.histGenotypeEnvelopeMaximum(maximum_of_genotype_hist_size,
                    maximum_of_genotype_hist_range);
            A_Ecosystem.histResourceAllocatedInCells(resource_alloc_hist_size,
                    resource_alloc_hist_range, resource_alloc_hist_start);
            A_Ecosystem.histResourceUptakenByCells(resource_uptake_hist_size,
                    resource_uptake_hist_range);
            A_Ecosystem.histsAgeAndResourceAtReprod(at_repr_intake_hist_size,
                    at_repr_intake_bin_range, AAR_hist_size, AAR_hist_range,
                    Reproduce_hist_size, Reproduce_hist_range, genome_size_histogram_size);
            A_Ecosystem.averageGenomeShape(environment_resolution);
        }
        A_Ecosystem.reproduceOrKillCell(mode_of_feeding, env, duplication_probability,
                deletion_probability, mutate_threshold, minimum_genotype_size,
                not_super_gene_param, gene_space_width, environment_resolution,
                one_over_total_env_space_surf, random_death_rate, minimal_cells_resource,
                reproduction_resourse_size, disparity_at_reproduction);
        A_Ecosystem.horizontalGeneTransfer(hgt_donor_prob, hgt_gene_trans);
    }
    A_Ecosystem.clearMutationCounters(half_numb_of_iterations);

    for (int k = half_numb_of_iterations + 1;  k <= number_of_iterations; k++) {
        // Tossing the next value of the environmental conditions 
        env = envGenerate(turbulence, twice_turbulence, p_rngEngine->rngMTgetDauble(), env);
        // all the stuff related to living
        A_Ecosystem.feedCells(mode_of_feeding, env, gene_replication_cost, core_genes,
                metabolic_costs, gain_from_the_genotype_ratio);
        // data extraction from the simulation
        if (k % sampling_interv == 0) {
            A_Ecosystem.writeGeneralDataToFile(k, env); // better not switch off this one
            A_Ecosystem.histGenomeSize(genome_size_histogram_size);
            A_Ecosystem.histCellsAge(AgeHistSize);
            A_Ecosystem.histGenotypeEnvelopeSurface(envelope_of_genotype_hist_size,
                    envelope_of_genotype_hist_range);
            A_Ecosystem.histGenotypeEnvelopeMaximum(maximum_of_genotype_hist_size,
                    maximum_of_genotype_hist_range);
            A_Ecosystem.histResourceAllocatedInCells(resource_alloc_hist_size,
                    resource_alloc_hist_range, resource_alloc_hist_start);
            A_Ecosystem.histResourceUptakenByCells(resource_uptake_hist_size,
                    resource_uptake_hist_range);
            A_Ecosystem.histsAgeAndResourceAtReprod(at_repr_intake_hist_size,
                    at_repr_intake_bin_range, AAR_hist_size, AAR_hist_range,
                    Reproduce_hist_size, Reproduce_hist_range, genome_size_histogram_size);
            A_Ecosystem.averageGenomeShape(environment_resolution);
        }
        A_Ecosystem.reproduceOrKillCell(mode_of_feeding, env, duplication_probability,
                deletion_probability, mutate_threshold, minimum_genotype_size,
                not_super_gene_param, gene_space_width, environment_resolution,
                one_over_total_env_space_surf, random_death_rate, minimal_cells_resource,
                reproduction_resourse_size, disparity_at_reproduction);
        A_Ecosystem.horizontalGeneTransfer(hgt_donor_prob, hgt_gene_trans);
    }
    A_Ecosystem.cumulativeMutationNumbers(number_of_iterations);
    A_Ecosystem.writeGenotypesToFile(env);
    AllTheGenomes.open("GenomesOfPopulation.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    AllTheGenomes << "# " << env << std::endl;
    //------The end of the main stuff------

    // ------ Post-production -------------
    p_rngEngine->release();
    pTagging_system->release();

    CellsAgeFile.close();
    time(&end);
    dif = difftime(end, start);


    CellsAgeFile.open("CellsAgeData.dat", std::ios::out | std::ios::ate | std::ios::app);
    CellsAgeFile << "# Running this program took " << dif << " seconds." << std::endl;
    CellsAgeFile.close();

    GetotypeHistFile.open("GenomeSizeData.dat", std::ios::out | std::ios::ate | std::ios::app);
    GetotypeHistFile << "# Running this program took " << dif << " seconds." << std::endl;
    GetotypeHistFile.close();

    GeneralDataFile.open("GeneralData.dat", std::ios::out | std::ios::ate | std::ios::app);
    GeneralDataFile << "# Running this program took " << dif << " seconds." << std::endl;
    GeneralDataFile.close();

    FrameSizeFile.open("FrameSizeData.dat", std::ios::out | std::ios::ate | std::ios::app);
    FrameSizeFile << "# Running this program took " << dif << " seconds." << std::endl;
    FrameSizeFile.close();

    FrameMaxFile.open("FrameMaxData.dat", std::ios::out | std::ios::ate | std::ios::app);
    FrameMaxFile << "# Running this program took " << dif << " seconds." << std::endl;
    FrameMaxFile.close();

    ResourceInCellsFile.open("RecourceInCells.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    ResourceInCellsFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    ResourceInCellsFile.close();

    ResourceUptakeCellsFile.open("RecourceUptakenCells.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    ResourceUptakeCellsFile << "# Running this program took " << dif << " seconds." \
      << std::endl;
    ResourceUptakeCellsFile.close();

    AgeAtFirstReproductionFile.open("AgeTillReproduction.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    AgeAtFirstReproductionFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    AgeAtFirstReproductionFile.close();

    NumberOfReproductionsFile.open("NumberOfReproductions.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    NumberOfReproductionsFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    NumberOfReproductionsFile.close();

    AverageIntakeTillReprodFile.open("AvarIntakeAtRepr.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    AverageIntakeTillReprodFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    AverageIntakeTillReprodFile.close();


    AverageTotGainTillReprodFile.open("AvarGainAtRepr.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    AverageTotGainTillReprodFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    AverageTotGainTillReprodFile.close();

    GeneNumberOfReprodFile.open("GeneNumberOfRepr.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    GeneNumberOfReprodFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    GeneNumberOfReprodFile.close();

    MutationRecords.open("MutationRecords.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    MutationRecords << "# Running this program took " << dif << \
        " seconds." << std::endl;
    MutationRecords.close();
    
    AllTheGenomes.open("GenomesOfPopulation.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    AllTheGenomes << "# Running this program took " << dif << \
        " seconds." << std::endl;
    AllTheGenomes.close();
    
    histGenotypeFile.open("AvaregeGenotype.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    histGenotypeFile << "# Running this program took " << dif << \
        " seconds." << std::endl;
    histGenotypeFile.close();
    
    GenotEnvelMean.open("GenotypeEnvelMean.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    GenotEnvelMean << "# Running this program took " << dif << \
        " seconds." << std::endl;
    GenotEnvelMean.close();
    
    GenotEnvelSTD.open("GenotypeEnvelSTD.dat",
            std::ios::out | std::ios::ate | std::ios::app);
    GenotEnvelSTD << "# Running this program took " << dif << \
        " seconds." << std::endl;
    GenotEnvelSTD.close();
    
    return 0;
}
