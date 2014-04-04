#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Iterates through directories with different model runs results and collects
data about environmental conditions stored in GeneralData.dat and some of the
parameters from ModelParams.dat. Then compares difference between variation of
the environment within cells expected life span (calculated from the actual
death rate) and the overall mean for the whole run.


Created on Tue Aug 16 00:03:32 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import os, re
import linecache as ln

LabelFontSize=22
TickSize=19

def LoadEnvCond(arg, dirname, files):
       for file in files:
            Grand_mean = p.nan
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname,'GeneralData.dat'):
                data = p.genfromtxt(filepath)
                if data[-1,4] != 0.0:
                    grand_mean = data[:, 1].mean()
                    frame_scaling = 1.0 / (data[1, 0] - data[0, 0])
                    filepath_turb = os.path.join(dirname,'ModelParams.dat')
                    l = re.split(" ", ln.getline(filepath_turb, 6))
                    turb_param = float(l[6])
                    l = re.split(" ", ln.getline(filepath_turb, 29))
                    rand_death_factor = float(l[6])
                    frame_half_size = p.rint((1.0 / (data[:, 11] \
                        / data[:, 4]).mean() ) * 0.5 * frame_scaling)
                    mean_turbul = p.zeros((data.shape[0]-(2.0 \
                        * frame_half_size), 2))
                    mean_turbul[:,0] = data[frame_half_size : \
                        - frame_half_size][:, 0].copy()
                    k = frame_half_size
                    for j in range(mean_turbul.shape[0]):
                        mean_turbul[j, 1] = p.absolute(data[k \
                            - frame_half_size : k + frame_half_size, 1].mean() \
                            - grand_mean)
                        k = k + 1
                    Grand_mean = mean_turbul[:, 1].mean()
                    arg.append((Grand_mean, turb_param, rand_death_factor))
                else:
                    break

GrandMeansData = []
os.path.walk(os.getcwd(), LoadEnvCond, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[1])

TheMeans = p.zeros((len(GrandMeansData), 2))
i = 0
for item in GrandMeansData:
    TheMeans[i,0] = item[0]
    TheMeans[i,1] = item[1]
    i += 1
xx = '# %(num)1.6f'%{"num":GrandMeansData[0][2]}
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)] # removing NaN's from the output
p.savetxt("Predictability.dat", TheMeans, fmt='%1.6f')
doc = open('Predictability.dat', 'a')
doc.write(xx)

p.figure(9, figsize=(1000,600))
p.plot(TheMeans[:, 1], TheMeans[:, 0], 'ko-')
p.fill_between(TheMeans[:, 1],TheMeans[:, 0], 0.0, color=(0.75, 0.75, 0.75, 0.75))
p.axis([0, 0.5, 0, 0.70001])
p.xlabel('Turbulence level ($ T $)', fontsize=LabelFontSize)
p.ylabel("Average difference from \n grand mean of env conditions",
         fontsize=LabelFontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)
p.show()
