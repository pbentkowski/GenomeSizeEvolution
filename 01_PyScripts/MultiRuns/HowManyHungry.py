#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 00:08:48 2010

@author: piotr bentkowski
"""

import pylab as p
import os
import re
import linecache as ln

important_line = 22  # change here to pick an interesting parameter

FontSize = 22
TickSize = 18


def CountHunry(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'RecourceUptakenCells.dat'):
            uptake_res_data = p.genfromtxt(filepath)
            uptake_res_data = uptake_res_data[-p.floor(uptake_res_data.shape[0]):]
            for k in xrange(uptake_res_data.shape[0]):
                uptake_res_data[k, :] = uptake_res_data[k, :] \
                    / uptake_res_data[k, :].sum()
            means_uptake = uptake_res_data[:, 0].mean()
            STD_uptake = uptake_res_data[:, 0].std()
            filepath_turb = os.path.join(dirname, 'ModelParams.dat')
            l = re.split(" ", ln.getline(filepath_turb, 6))
            turb_param = float(l[6])
            l = re.split(" ", ln.getline(filepath_turb, important_line))
            interesting_param = float(l[6])
            arg.append((means_uptake, STD_uptake,
                        turb_param, interesting_param))

Hungry = []
os.path.walk(os.getcwd(), CountHunry, Hungry)
Hungry = sorted(Hungry, key=lambda data_sort: data_sort[2])
#print Hungry

#  ---- loading and cheking data ---
TheHungryOnes = p.zeros((len(Hungry), 4))
i = 0
for item in Hungry:
    TheHungryOnes[i, 0] = 1.0 - item[0]
    TheHungryOnes[i, 1] = item[1]
    TheHungryOnes[i, 2] = item[2]
    TheHungryOnes[i, 3] = item[3]
    i += 1

print TheHungryOnes
 # removing NaN's from the output
TheHungryOnes = TheHungryOnes[~p.isnan(TheHungryOnes).any(1)]

p.savetxt("HugryCells.dat", TheHungryOnes)

p.figure(10, figsize=(12, 6))
#p.subplot(221)
p.plot(TheHungryOnes[:, 2], TheHungryOnes[:, 0], 'ko-')
p.fill_between(TheHungryOnes[:, 2], TheHungryOnes[:, 0] + TheHungryOnes[:, 1],
               TheHungryOnes[:, 0] - TheHungryOnes[:, 1],
               color=(0.75, 0.75, 0.75, 0.75))
p.axis([0, 0.5, 0, 0.8001])
#p.title(r'$\delta = %(delt)1.3f $' % {"delt": TheHungryOnes[0, 3]},
#        fontsize=FontSize)
p.xlabel('turbulence level ($ T $)', fontsize=FontSize)
p.ylabel('fraction of cells that were fed', fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)
p.show()
