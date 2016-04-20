#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Opens the file MutationStats.dat to an array and just plots what's in it. Use
script Mutation_rates.py to create the MutationStats.dat file from the model
runs results.
This is HGT version (modified on 03/04/2012)

Created on Thu Jul 14 17:29:27 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import Y_axis_div as div  # this one is mine

AxisLabelFontSize = 22
AxisTickFontSize = 22
AnnotateFontSize = 19

SectOnYaxis = 5
SectOnYaxis_2 = 5

TheFinalArray = p.genfromtxt("MutationStats.dat")

max_modif = (TheFinalArray[:, 3] + TheFinalArray[:, 4]).max()
max_dupl = (TheFinalArray[:, 5] + TheFinalArray[:, 6]).max()
max_del = (TheFinalArray[:, 7] + TheFinalArray[:, 8]).max()
max_of_maxes = max(max_modif, max_dupl, max_del)
[MIN, MaxMutationNumber, interv] = div.divide_Y_axis(0.0, max_of_maxes,
                                                     SectOnYaxis)

Letter_location = MaxMutationNumber - interv / 2.0

p.figure(1, figsize=(1000, 600))
tick_marks = p.arange(MIN, MaxMutationNumber + 0.001, interv)
p.subplot(221)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 3], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 3] - TheFinalArray[:, 4],
               TheFinalArray[:, 3] + TheFinalArray[:, 4],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'A', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of gene modifications', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim(xmin=0.0)
p.grid(True)

p.subplot(222)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 5], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 5] - TheFinalArray[:, 6],
               TheFinalArray[:, 5] + TheFinalArray[:, 6],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'B', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of duplications', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim(xmin=0.0)
p.grid(True)

p.subplot(223)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 7], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 7] - TheFinalArray[:, 8],
               TheFinalArray[:, 7] + TheFinalArray[:, 8],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'C', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of deletions', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim(xmin=0.0)
p.grid(True)

max_clone_num = max(TheFinalArray[:, 2])
[MIN, MaxClonalStrain, intervvv] = div.divide_Y_axis(0.0, max_clone_num,
                                                     SectOnYaxis_2)

p.subplot(224)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 9], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 9] - TheFinalArray[:, 10],
               TheFinalArray[:, 9] + TheFinalArray[:, 10],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'D', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of HGT events', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim(xmin=0.0)
p.grid(True)

max_clone_num = max(TheFinalArray[:, 2])
[MIN, MaxClonalStrain, intervvv] = div.divide_Y_axis(0.0, max_clone_num,
                                                     SectOnYaxis_2)

p.figure(2, figsize=(1000, 600))
tick_marks_04 = p.arange(0., MaxClonalStrain + 0.1, intervvv)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 2], 'ko-')
ax = p.text(0.475, intervvv / 2.0, 'E', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of clonal strains', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_04, size=AxisTickFontSize)
p.grid(True)

p.show()
