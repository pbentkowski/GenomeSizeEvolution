#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:11:19 2010
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import os
import re
import linecache as ln

AxisLabelFontSize = 15
AxisTickFontSize = 15
AnnotateFontSize = 14


def LoadEnvelData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'Evelopes.dat'):
            envelopes = p.genfromtxt(filepath)
            # goes like this: mean envelope surface, turbulance level T, STD of
            # envelope surface
            arg.append((envelopes[:, 0], envelopes[:, 2], envelopes[:, 1]))


def LoadGeneNumbers(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GenomeSizes.dat'):
            genes = p.genfromtxt(filepath)
            # goes like this: mean number of genes, turbulance level T, STD of
            # number of genes
            arg.append((genes[:, 0], genes[:, 2], genes[:, 1]))


def LoadIterestParam(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'Evelopes.dat'):
                l = re.split(" ", ln.getline(filepath, 1))
                gene_cost = float(l[1])
                arg.append((gene_cost))

EnvelData = []
os.path.walk(os.getcwd(), LoadEnvelData, EnvelData)
GeneData = []
os.path.walk(os.getcwd(), LoadGeneNumbers, GeneData)
Tags = []
os.path.walk(os.getcwd(), LoadIterestParam, Tags)

p.figure(1, figsize=(10, 6))
p.subplot(221)
i = 0
for item in EnvelData:
    XX = float(item[1][-i-1])
    YY = float(item[0][-3])
    ax = p.plot(item[1], item[0], 'k-')
    ax = p.annotate('%1.3f' % (Tags[i], ), xy=(XX, YY),  xycoords='data',
                    fontsize=AnnotateFontSize)
    i = i + 1
#p.xlabel('turbulence level ($ T $)', fontsize = AxisLabelFontSize)
p.ylabel('mean genotype\'s envelope surface', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()

p.subplot(223)
i = 0
for item in GeneData:
    XX = float(item[1][-i-1])
    YY = float(item[0][-3])
    ax = p.plot(item[1], item[0], 'k-')
    ax = p.annotate('%1.3f' % (Tags[i],), xy=(XX, YY),  xycoords='data',
                    fontsize=AnnotateFontSize)
    i = i + 1
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.ylabel('mean number of genes', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()

p.subplot(122)
for j in xrange(len(EnvelData)):
    for i in xrange(len(EnvelData[j][0])):
#        ax = p.errorbar(GeneData[j][0][i], EnvelData[j][0][i],
#                         xerr=GeneData[j][2][i], yerr=EnvelData[j][2][i],
#                         fmt='o', ecolor='k')
        ax = p.plot(GeneData[j][0][i], EnvelData[j][0][i], 'ko')
p.xlabel('number of genes', fontsize=AxisLabelFontSize)
p.ylabel('surface of genotype\'s envelope to env. ratio',
         fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()
p.show()
