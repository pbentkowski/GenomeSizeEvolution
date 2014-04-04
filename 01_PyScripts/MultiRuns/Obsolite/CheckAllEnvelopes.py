#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 21:41:03 2010

@author: piotr
"""

import os
import re
import pylab as p
import linecache as ln


def LoadMyData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'Evelopes.dat'):
            envelopes = p.genfromtxt(filepath)
            arg.append((envelopes[:, 0], envelopes[:, 2]))


def LoadGeneCost(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'Evelopes.dat'):
                l = re.split(" ", ln.getline(filepath, 1))
                gene_cost = float(l[1])
                arg.append((gene_cost))

AllData = []
os.path.walk(os.getcwd(), LoadMyData, AllData)
Tags = []
os.path.walk(os.getcwd(), LoadGeneCost, Tags)

#print AllData[1][1][1]

p.figure(1, figsize=(1000, 600))
i = 0
for item in AllData:
    XX = float(item[1][-i-3])
    YY = float(item[0][33])
    ax = p.plot(item[1], item[0], 'k-')
    ax = p.annotate('%4.2f' % (Tags[i],), xy=(XX, YY), xycoords='data',
                    sfontsize=14)
    i = i + 1
p.xlabel('turbulence level ($ T $)', fontsize=17)
p.ylabel('mean genotype\'s envelope surface', fontsize=17)
p.xticks(size=16)
p.yticks(size=16)
p.grid()
p.show()
