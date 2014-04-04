#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Iterates through directories with data and calculates the mean death rate for
stable part of the models runs. Puts then in a figure.

Created on Fri Feb 17 17:55:32 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import os, re
import linecache as ln

important_line = 53 # change here to pick an interesting parameter

FontSize = 22
TickSize = 18

def CountTrueDeathRate(arg, dirname, files):
       for file in files:
            Grand_mean = p.nan
            Grand_STD = p.nan
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'GeneralData.dat'):
                data = p.genfromtxt(filepath)
                if data[-1, 4] != 0.0:
                    half_of_data = p.floor(data.shape[0] / 2)
                    data_chopped = data[half_of_data:-1, :]
                    Grand_mean = (data_chopped[:, 11] / data_chopped[:, 4]).mean()
                    Grand_STD = (data_chopped[:, 11]/ data_chopped[:, 4]).std()
                    filepath_turb = os.path.join(dirname, 'ModelParams.dat')
                    l = re.split(" ", ln.getline(filepath_turb, 6))
                    turb_param = float(l[6])
                    l = re.split(" ", ln.getline(filepath_turb, important_line))
                    interesting_param = float(l[6])
                    arg.append((Grand_mean, Grand_STD, turb_param,
                                interesting_param))
                else:
                    break

GrandMeansData = []
os.path.walk(os.getcwd(), CountTrueDeathRate, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])
#  ---- loading and checking data ---
TheMeans = p.zeros((len(GrandMeansData), 4))
i = 0
for item in GrandMeansData:
    TheMeans[i,0] = item[0]
    TheMeans[i,1] = item[1]
    TheMeans[i,2] = item[2]
    TheMeans[i,3] = item[3]
    i += 1
xx = '# %(num)1.4f'%{"num":GrandMeansData[0][3]}

print TheMeans
print xx
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)] # removing NaN's from the output

p.savetxt("DeathRates.dat", TheMeans, fmt='%1.5f')
doc = open('DeathRates.dat', 'a')
doc.write(xx)

p.figure(1, figsize=(1000,600))
#p.subplot(221)
p.plot(TheMeans[:,2],TheMeans[:,0], 'ko-' )
p.fill_between(TheMeans[:,2], TheMeans[:,0]+TheMeans[:,1],
               TheMeans[:,0]-TheMeans[:,1], color=(0.75, 0.75, 0.75, 0.75))
p.axis([0, 0.5, 0, 0.014])
#p.title(r'$\delta = %(delt)1.3f $'% {"delt": TheMeans[0,3]}, fontsize=FontSize)
p.xlabel('turbulence level ($ T $)', fontsize=FontSize)
p.ylabel('death rate', fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)
p.show()