#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots the genomes' sizes and the ratio of the genomes' envelopes to the surface
of the whole environment for number of runs with different turbulence levels.
Uses files Evelopes.dat and GenomeSizes.dat computed by scripts
MeansEnvelope.py and MeansGenomeSize.py respectively.

Created on Thu Oct 7 14:42:53 2010
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import sys

x_range = 0.5  # range of the x-axis in something vs turbulance plots
FontSize_axis = 20
AxisTickFontSize = 17
if len(sys.argv) > 1:
    titl = str(sys.argv[1])
else:
    titl = ''

envelopes = p.genfromtxt("Evelopes.dat")
genome_sizes = p.genfromtxt("GenomeSizes.dat")

if p.size(envelopes) != p.size(genome_sizes):
    print "Number of entries in Evelopes.dat and GenomeSizes.dat do not",
    print " match! Check those files"

p.figure(9, figsize=(18, 12))
p.subplot(322)
p.title(titl, fontsize=24)
y_range_01 = 20.01
tick_marks_01 = p.arange(0., y_range_01, 4.)
p.plot(genome_sizes[:, 2], genome_sizes[:, 0], 'ko-')
p.fill_between(genome_sizes[:, 2], genome_sizes[:, 0] + genome_sizes[:, 1],
               genome_sizes[:, 0] - genome_sizes[:, 1],
               color=(0.75, 0.75, 0.75, 0.75))
p.xlabel('turbulence level ($ T $)', fontsize=FontSize_axis)
p.ylabel('number of genes', fontsize=FontSize_axis)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_01, size=AxisTickFontSize)
p.axis([0, x_range, 0, y_range_01])
p.grid(True)

p.subplot(324)
y_range_02 = 0.6001
tick_marks_02 = p.arange(0., y_range_02, 0.1)
p.plot(envelopes[:, 2], envelopes[:, 0], 'ko-')
p.fill_between(envelopes[:, 2], envelopes[:, 0] + envelopes[:, 1],
               envelopes[:, 0] - envelopes[:, 1],
               color=(0.75, 0.75, 0.75, 0.75))
p.xlabel('turbulence level ($ T $)', fontsize=FontSize_axis)
p.ylabel('genotype\'s \n envelope size', fontsize=FontSize_axis)
p.axis([0, x_range, 0, y_range_02])
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_02, size=AxisTickFontSize)
p.grid(True)

p.subplot(326)
p.errorbar(genome_sizes[:, 0], envelopes[:, 0], xerr=genome_sizes[:, 1],
           yerr=envelopes[:, 1], fmt='o', ecolor='k')
p.plot(genome_sizes[:, 0], envelopes[:, 0], 'ko')
p.xlabel('number of genes', fontsize=FontSize_axis)
p.ylabel('genotype\'s  \n envelope size', fontsize=FontSize_axis)
p.xticks(tick_marks_01, size=AxisTickFontSize)
p.yticks(tick_marks_02, size=AxisTickFontSize)
p.axis([0, y_range_01, 0, y_range_02])
p.grid(True)

p.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
p.show()
