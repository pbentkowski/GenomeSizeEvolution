#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Put a doc string here...

Created on Mon May 23 18:37:03 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import os
import re
import pylab as p
import linecache as ln
# the one below is mine
import file_len as fl

AxisLabelFontSize = 22
AxisTickFontSize = 22
AnnotateFontSize = 19

gene_num_lim = 40
gene_surf_lim = 0.7000001
dec_places = '%1.4f'


def LoadEnvelData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'Evelopes.dat'):
            envelopes = p.genfromtxt(filepath)
            # goes like this: mean envelope surface, turbulance level T,
            # STD of envelope surface
            arg.append((envelopes[:, 0], envelopes[:, 2], envelopes[:, 1]))


def LoadGeneNumbers(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GenomeSizes.dat'):
            genes = p.genfromtxt(filepath)
            # goes like this: mean number of genes, turbulance level T,
            # STD of number of genes
            arg.append((genes[:, 0], genes[:, 2], genes[:, 1]))


def LoadGeneCost(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'Evelopes.dat'):
                l = re.split(" ", ln.getline(filepath, fl.file_len(filepath)))
                death_rate = float(l[1])
                arg.append((death_rate))

EnvelData = []
os.path.walk(os.getcwd(), LoadEnvelData, EnvelData)
GeneData = []
os.path.walk(os.getcwd(), LoadGeneNumbers, GeneData)
Tags = []
os.path.walk(os.getcwd(), LoadGeneCost, Tags)


p.figure(1, figsize=(15, 9))
p.subplot(221)
i = 0
for item in EnvelData:
    XX = 0.05 + i/15.0  # float(item[1][-(i*4)-3])
    YY = float(item[0][-3])
    ax = p.plot(item[1], item[0], 'k-')
    ax = p.annotate(dec_places % (Tags[i],), xy=(XX, YY),  xycoords='data',
                    fontsize=AnnotateFontSize)
    i = i + 1
ax = p.text(0.45, 0.05, 'A', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=24)
#p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.ylabel('mean genotype\'s \n envelope surface', fontsize=AxisLabelFontSize)
p.axis([0, 0.5, 0, gene_surf_lim])
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()

p.subplot(223)
i = 0
for item in GeneData:
    XX = 0.05 + i/15.0  # float(item[1][10 + i])
    YY = float(item[0][-3])
    ax = p.plot(item[1], item[0], 'k-')
    ax = p.annotate(dec_places % (Tags[i],), xy=(XX, YY),  xycoords='data',
                    fontsize=AnnotateFontSize)
    i = i + 1
ax = p.text(0.45, 2.5, 'B', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=24)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.ylabel('mean number of genes', fontsize=AxisLabelFontSize)
p.axis([0, 0.5, 0, gene_num_lim])
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()

p.subplot(122)
for j in xrange(len(EnvelData)):
    for i in xrange(len(EnvelData[j][0])):
#        ax = p.errorbar(GeneData[j][0][i], EnvelData[j][0][i],
#                         xerr=GeneData[j][2][i], yerr=EnvelData[j][2][i],
#                        fmt='o', ecolor='k')
        ax = p.plot(GeneData[j][0][i], EnvelData[j][0][i], 'ko')
ax = p.text(gene_num_lim - 2.5, 0.05, 'C', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=24)
p.xlabel('number of genes', fontsize=AxisLabelFontSize)
p.ylabel('surface of genotype\'s envelope to env. ratio',
         fontsize=AxisLabelFontSize)
p.axis([0, gene_num_lim, 0, gene_surf_lim])
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()

p.show()
