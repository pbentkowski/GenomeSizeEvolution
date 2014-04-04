//      gene.cpp
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
#include <sstream>
#include <iomanip>
#include <gsl/gsl_math.h>
#include "rngEngine.h"

/** \brief Core method. Constructor - sets a gene, but with no values. 
 * 
 * Just a constructor.
 * 		 		 		 		 		 		   
 */
Gene::Gene() {
}

/** \brief Core method. Sets values of a gene.
 * 
 * Assigns the values of genes by tossing each value. Keeps the surface
 * under the Gaussian curve fixed.
 * 
 * @param not_super_gene_param - surface under the Gaussian curve
 * @param gene_space_width - constant helping to calculate growth rate
 *
 */
void Gene::setGeneValues(double not_super_gene_param, double gene_space_width) {
    // GeneValue[0] - hight ; GeneValue[1] - variance ; GeneValue[2] - mean
    // of the Gaussian function
    rngEngine* p_rngEngine = rngEngine::getInstance();
    // hight [0,1)
    GeneValue[0] = p_rngEngine->rngMTgetDauble();
    // variance; it scales to hight to lover the probability of emergence
    // of a 'supergene'
    GeneValue[1] = not_super_gene_param / (GeneValue[0] * sqrt(2.0 * M_PI));
    // mean
    GeneValue[2] = 2.0 * gene_space_width * p_rngEngine->rngMTgetDauble() \
        - gene_space_width;
}

/** \brief Core method. Sets new values of genes if tossed below mutation threshold.  
 * 
 * Assigns the values of genes by tossing  if random factor is below fixed
 * mutation threshold.
 * @param not_super_gene_param -surface under the Gaussian curve
 * @param mutate_threshold - gene modification probability (mutation)
 * @param gene_space_width - the width of the env space genes can take
 */
void Gene::mutate(double not_super_gene_param, double mutate_threshold,
        double gene_space_width) {
    rngEngine* p_rngEngine = rngEngine::getInstance();
    // mutating hight
    if (p_rngEngine->rngMTgetDauble() < mutate_threshold) {
        // hight, [0, 1)
        GeneValue[0] = p_rngEngine->rngMTgetDauble();
    }
    // mutating variance
    if (p_rngEngine->rngMTgetDauble() < mutate_threshold) {
        // variance; it scales to hight to lover the probability of emergence
        // of a 'supergene'
        GeneValue[1] = not_super_gene_param / (GeneValue[0] * sqrt(2.0 * M_PI));
    }
    // mutating mean
    if (p_rngEngine->rngMTgetDauble() < mutate_threshold) {
        // mean
        GeneValue[2] = 2.0 * gene_space_width * p_rngEngine->rngMTgetDauble() \
            - gene_space_width;
    }
}

/** \brief Core method. Calculates gross growth rate from gene.
 * 
 * Calculates  \f$ u(x,c,\sigma,A) = A  \exp( \frac{-(x-c)^2}{2\sigma^2} )\f$
 * for value of \f$x\f$ given as argument. Uses
 * <a href="http://www.gnu.org/software/gsl/manual/html_node/Small-integer-powers.html">
 * GNU Scientific Library</a>
 * funcion to calculate power of 2.
 * @param Env - current value of the environmental conditions.
 * @return - growth rate given by the gene for this env condition (double)
 */
double Gene::calculateGrowthFromGene(double Env) {
    return (GeneValue[0] * exp((-gsl_pow_2(Env - GeneValue[2])) \
        / (2 * gsl_pow_2(GeneValue[1]))));
}

/** \brief Data collecting method. Converts all three gene's values into a
 * single sting. 
 * 
 * By changing the number in std::setprecision(n) you can change number of displayed
 * digits. May be important as it influences the length of a line in some data 
 * output files.
 * 
 * @param void
 * @return three values of the gene in a form of STD string 
 * <a href="http://www.cplusplus.com/reference/string/string/">STD string</a>.
 */
std::string Gene::geneToString() {
    std::ostringstream gene; 
    gene << std::setprecision(4) << GeneValue[0] << " " << GeneValue[1] << " " \
            << GeneValue[2];
    return gene.str();
}

/** \brief Data collecting method. Returns just one chosen value of a gene.
 * 
 *  
 * Takes index number of gene array and returns value of corresponding gene
 * parameter. If index is out of gene array range prints an error message and
 * returns -1.
 * @param GeneValueIndex - index of one of the Gaussian curve params (can be 0, 1 or 2)
 * @return GeneValue[GeneValueIndex] - value of the chosen index (double) or -1 (int)
 */
double Gene::getOneValue(unsigned int GeneValueIndex) {
    if ((sizeof (GeneValue) / sizeof (GeneValue[0])) <= GeneValueIndex) {
        std::cout << "GeneValueIndex provided in getOneValue(int) is \
            outside GeneValues array range" << std::endl;
        return -1;
    } else {
        return GeneValue[GeneValueIndex];
    }
}
