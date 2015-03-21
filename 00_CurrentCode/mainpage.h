/** 
 * @mainpage "Modelling Genome Size Constrains"
 * 
 * @author <b>Piotr Bentkowski</b>: \n <a href="mailto:bentkowski.piotr@gmail.com
 * ?subject=Genome Streamlining Model">bentkowski.piotr@gmail.com</a> \n
 * \n
 * Developed during Piotr's studies at:\n
 * Laboratory for Global Marine and Atmospheric Chemistry \n
 * School of Environmental Sciences \n
 * University of East Anglia \n
 * NR4 7TJ, Norwich \n
 * United Kingdom \n
 * \n
 * Under supervision of Dr. Thomas Mock \n
 * With help from: Dr. Cock van Oosterhout, Dr. Hywel Williams, Prof. Tim Lenton \n
 * @version 1.1
 * 
 * @section intro Introduction
 * This program is an implementation of a model of the evolution of genome size described
 * in Piotr Bentkowski's Ph.D. thesis:\n
 * <a href="https://ueaeprints.uea.ac.uk/50553/">Modelling evolution of genome size in prokaryotes 
 * in response to changes in their abiotic environment</a> (2014) defended at the University of East Anglia. \n
 * 
 * @section Licence
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version. \n
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. \n
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * @section Compilation
 * The program has some dependencies on
 * <a href="http://www.gnu.org/software/gsl/"> GNU Scientific Library</a>. Should compile
 * smoothly on most GNU/Linux distros with GSL installed.
 *  Compilation works fine on Ubuntu 12.04 LTS and 14:04 LTS with a command: \n 
 * \n 
 *<b>
 * g++ -O2 -o TheModel -L/usr/local/lib TheModel.cpp genotype.cpp gene.cpp
 * cell.cpp ecosystem.cpp rngEngine.cpp tagging_system.cpp -lgsl -lgslcblas -lm
 * </b>\n 
 * \n
 * It also worked fine on UEA's GRACE Computational Cluster compiled with Intel C++ Compiler (icc) 
 * 
 * 
 * @section Parameters
 * <b>Program takes exactly 45 parameters. These are:</b> \n
 * <b> 00 </b>- program's name \n
 * <b> 01 </b>- the turbulence level (the length of the step in the random walk; 
 * between 0 and 1) \n
 * <b> 02 </b>- resolution of the environment (how dense the environment space is
 *  being sampled, well < 1) \n
 * <b> 03 </b>- the value = 1/[env space surface] (calculate it on a sheet of paper 
 * cos it's gonna be faster that way)  \n
 * <b> 04 </b>- a seed for the Mersenne Twister random number generator (a positive
 *  integer) \n
 * <b> 05 </b>- the number of the iterations in the model's run (a positive
 *  integer) \n
 * <b> 06 </b>- the 'not-super-gene' parameter (surface under the Gaussian curve)\n
 * <b> 07 </b>- the probability of mutation of existing gene (between 0 and 1) \n
 * <b> 08 </b>- the fraction of env space is occupied by the gene space (suggested 
 * value = 1) \n
 * <b> 09 </b>- the probability of deleting a gene (between 0 and 1) \n
 * <b> 10 </b>- the probability of duplicating a gene (between 0 and 1) \n
 * <b> 11 </b>- the maximal number of genes at model's initialisation (don't aim 
 * to high, rather dozens then hundreds) \n
 * <b> 12 </b>- the minimal number of genes at initialisation (must be lower then
 *  #11)\n
 * <b> 13 </b>- the minimal number of genes at all times (at least 1 must be 
 * present or it may result in a segmentation fault) \n
 * <b> 14 </b>- the number of core genes (might be 0)\n
 * <b> 15 </b>- the cost of maintaining one gene (>= 0)\n
 * <b> 16 </b>- the gain from maximal fitness (>= 0 , maximal amount of the 
 * resource a cell can gain in one go) \n
 * <b> 17 </b>- the constant metabolic cost (>= 0, a cost of being a cell and
 *  being alive)\n
 * <b> 18 </b>- the minimal resource below which a cell dies\n
 * <b> 19 </b>- the minimal resource above which a cell divides\n
 * <b> 20 </b>- the maximal resource bonus given to a cell at the initialisation 
 * of the model (actual value is being randomly selected between 0 and this 
 * number) \n
 * <b> 21 </b>- the probability of a random death (between 0 and 1)\n
 * <b> 22 </b>- the total amount of resource the ecosystem has (a pretty big number,
 *  please)\n
 * <b> 23 </b>- the number of bins in the genome size histogram \n
 * <b> 24 </b>- the number of bins in the histogram of the envelope of 
 * genome size \n
 * <b> 25 </b>- the width of the bin in the histogram of the envelope of 
 * genome size \n
 * <b> 26 </b>- the number of bins in the histogram of maximums of the fitness
 * function \n
 * <b> 27 </b>- the width of bin in the histogram maximums of the fitness
 * function \n
 * <b> 28 </b>- the number of bins in the histogram of the age of cells \n
 * <b> 29 </b>- the number of bins in the histogram of the resource allocated
 * to cells\n
 * <b> 30 </b>- the width of bins in the histogram of the resource allocated
 * to cells \n
 * <b> 31 </b>- the beginning value of the histogram of the resource allocated
 * to cells \n
 * <b> 32 </b>- the number of bins in the histogram of the resource up-taken
 * by cells in population \n
 * <b> 33 </b>- the width of bins in the histogram of the resource up-taken
 * by cells in population\n
 * <b> 34 </b>- the mode of feeding ('s' -sequential, 'p' -proportional; other
 * values will cause the program to prompt an error message and quit)\n
 * <b> 35 </b>- the parameter how uneven is the resource distribution when a cell
 * is dividing (has meaning only when program works with proportional feeding
 * mode; value between 0 and 0.5) \n
 * <b> 36 </b>- the number of bins in the histogram of the age at the first
 * reproduction \n
 * <b> 37 </b>- the width of bins in the histogram of the age at the first
 * reproduction \n
 * <b> 38 </b>- the number of bins in the histogram of the number of the reproduction
 * events per cell \n
 * <b> 39 </b>- the width of bins in the histogram of the number of the reproduction
 * events per cell \n
 * <b> 40 </b>- the number of bins of the histogram of the average intake per unit
 * of time until reproduction \n
 * <b> 41 </b>- the bin width of the histogram of the average intake per unit of
 * time until reproduction \n 
 * <b> 42 </b>- probability that a cell will become a gene donor in HGT \n
 * <b> 43 </b>- probability that a gene will get transferred from the donor's
 * genotype to recipient's genotype \n
 * <b> 44 </b>- time interval (steps) in which statistics are saved to data files\n
 * \n
 * This program can recognise simple errors in the argument list (e.g. value out of
 * range, negative values when positive are needed etc.), but will not recognise
 * when they don't make a 'biological' sense.
 * 
 * @section The Output
 * Program produces a number of text files containing desired data. The file
 * <i>ModelParams.dat</i> contains the values of parameters fed to the model. The
 * file <i>GenomeSizeData.dat</i> contains the main statistics of the model,
 * including the time stamp used in the visualisation process of the other files
 * with data. This two files are the basics to analyse simulation output. Other
 * files are rather optional and they can be switched off by commenting out some
 * lines in the main function (<i>TheModel.cpp</i> file) what should speed up
 * the program runtime.
 * \n\n
  *  Important files:\n
 * <b><i>ModelParams.dat</i></b> - contains run's parametrisation \n 
 * <b><i>GeneralData.dat</i></b> - contains basic statistics of the output \n
 * \n
 * Each line of these files represents one time step (in sync with <i>GeneralData.dat</i>): \n
 * <b><i>GenomeSizeData.dat</i></b> - number of genes in genomes; histograms \n
 * <b><i>GeneNumberOfRepr.dat</i></b> - number of genes in genomes but only of cells which reproduced; histograms \n
 * <b><i>CellsAgeData.dat</i></b> - age of the cells in population; histograms \n
 * <b><i>AgeTillReproduction.dat</i></b> - time passed since the last time cell reproduced (or was born); histograms \n
 * <b><i>FrameSizeData.dat</i></b> - size of the genotype's frame; histograms \n
 * <b><i>FrameMaxData.dat</i></b> - maximums of the genotype's frame; histograms \n 
 * <b><i>AvarGainAtRepr.dat</i></b> - cells' average total gain before reproduction; histograms \n
 * <b><i>AvarIntakeAtRepr.dat</i></b> - cells' average intake before reproduction; histograms \n
 * <b><i>RecourceInCells.dat</i></b> - resources allocated into cells in one time step; histograms \n
 * <b><i>NumberOfReproductions.dat</i></b> - how many times a cell reproduced; histograms \n
 * <b><i>RecourceUptakenCells.dat</i></b> - resources taken by cells in one time step; histograms \n
 * <b><i>AvaregeGenotype.dat</i></b> - U(x) values (intake values) in all the cells in the population; histograms \n
 * <b><i>AvarIntakeAtRepr.dat</i></b> - cells' average intake before reproduction; histograms \n
 * <b><i>FrameMaxData.dat</i></b> - maximums of the genotypes' frame; histograms \n
 * <b><i>GeneNumberOfRepr.dat</i></b> - cells' gene number distribution, only cells which reproduced; histograms \n
 * <b><i>GenotypeEnvelMean.dat</i></b> - mean values of U(x) function (the genome's envelope) \n
 * <b><i>GenotypeEnvelSTD.dat</i></b> - SD of U(x) function (the genome's envelope) \n
 * 
 * \n 
 * Other data files: \n
 * <b><i>MutationRecords.dat</i></b> - cumulative number of mutations accumulated by each living cell 
 * in a defined time period. There is a dedicated function Ecosystem::clearMutationCounters(int current_iteration) 
 * for zeroing number of mutation  \n
 * <b><i>GenomesOfPopulation.dat</i></b> - these are all the genomes in the ecosystem \n
 * 
 * @section Visualisation
 * Visualisation is done using Python scripts containing calls to Pylab library
 * (Linux distros will have it in its repos, for Windows you might check 
 * the <a href="http://code.google.com/p/pythonxy/">Python(x,y) project</a>). These scripts
 * can be found in <i>01_PyScripts</i> directory.
 */