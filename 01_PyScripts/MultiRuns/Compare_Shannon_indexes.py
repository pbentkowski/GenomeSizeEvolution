#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots a comparison of Shannon indexes saved in Shannons.dat files in different
directories belonging to different model runs. Uses a parameter determined in
the last line of Shannons.dat to label the plots.

Created on Thu Jul  7 15:38:59 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import os
import re
import pylab as p
import linecache as ln
# the one below is mine
import file_len as fl

AxisLabelFontSize = 23
AxisTickFontSize = 23
AnnotateFontSize = 20


def LoadShannonsData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'Shannons.dat'):
            shannons = p.genfromtxt(filepath)
            # goes like this: mean Shannon index, turbulance level T,
            # STD of Shannon index
            arg.append((shannons[:, 0], shannons[:, 2], shannons[:, 1]))


def LoadInterestinParam(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'Shannons.dat'):
                l = re.split(" ", ln.getline(filepath, fl.file_len(filepath)))
                fancy_param = float(l[1])
                arg.append((fancy_param))


ShannonData = []
os.path.walk(os.getcwd(), LoadShannonsData, ShannonData)
Tags = []
os.path.walk(os.getcwd(), LoadInterestinParam, Tags)

p.figure(1, figsize=(10, 6))
i = 0
for item in ShannonData:
    XX = float(item[1][-i-3])
    YY = float(item[0][-3])
    ax = p.plot(item[1], item[0], 'ko-', markersize=3)
    ax = p.fill_between(item[1], item[0] + item[2], item[0] - item[2],
                        color=(0.75, 0.75, 0.75, 0.75))
    ax = p.annotate('%1.3f' % (Tags[i],), xy=(XX, YY),  xycoords='data',
                    fontsize=AnnotateFontSize)
    i = i + 1
#p.xlabel('turbulence level ($ T $)', fontsize = AxisLabelFontSize)
p.ylabel('mean Shannon index ($ H $)', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()
p.show()
