#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 00:45:08 2010

@author: piotr
"""

import os
import pylab as p


def LoadMyEnvel(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'Evelopes.dat'):
            envelopes = p.genfromtxt(filepath)
            arg.append((envelopes[:, 0], envelopes[:, 1]))


def LoadMyGene(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GenomeSizes.dat'):
            genes = p.genfromtxt(filepath)
            arg.append((genes[:, 0], genes[:, 1]))


def main():
    """ """
    EnvelData = []
    os.path.walk(os.getcwd(), LoadMyEnvel, EnvelData)

    GeneData = []
    os.path.walk(os.getcwd(), LoadMyGene, GeneData)

    p.figure(11, figsize=(12, 6))
    for j in xrange(len(EnvelData)):
        for i in xrange(len(EnvelData[j][0])):
            ax = p.errorbar(GeneData[j][0][i], EnvelData[j][0][i],
                            xerr=GeneData[j][1][i], yerr=EnvelData[j][1][i],
                            fmt='o', ecolor='k')
    p.xlabel('number of genes', fontsize=17)
    p.ylabel('surface of genotype\'s envelope to env. ratio', fontsize=17)
    p.xticks(size=16)
    p.yticks(size=16)
    p.grid()
    p.show()

    print 'done'

if __name__ == "__main__":
    main()
