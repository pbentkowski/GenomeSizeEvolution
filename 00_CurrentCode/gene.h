//      gene.h
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


#include <iostream>

#ifndef _GENE_H_
#define _GENE_H_

/** 
 * @class Gene
 *
 * @brief Single gene class. Stores 3 double numbers necessary to draw
 * Gaussian curve in environmental space. 
 * Has methods to change them.
 */
class Gene {
public:
    // gene values
    /** 
     * @brief Array storing three parameters for plotting gene.
     * 
     * GeneValue[0] - hight ; GeneValue[1] - variance ; GeneValue[2] -
     * mean of Gaussian function in environmental space.
     */
    double GeneValue[3];
    // initialisation functions
    Gene();
    void setGeneValues(double not_super_gene_param, double gene_space_width);
    void mutate(double not_super_gene_param, double mutate_threshold,
            double gene_space_width);
    double calculateGrowthFromGene(double Env);
    std::string geneToString();
    double getOneValue(unsigned int GeneValueIndex);
};

#endif
