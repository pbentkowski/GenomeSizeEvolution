#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Iterates through the directories with data. Calculates the variance of the env
conditions for the whole time space and a mean variance for time frames equal
to the expected life span of the cell (calculated from the actual death rate)

Created on Sun Feb 26 00:47:21 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""

import pylab as p
import os, re
import linecache as ln

LabelFontSize = 20
TickSize = 17

def LoadEnvCond(arg, dirname, files):
       for file in files:
            Grand_var = p.nan
            Grand_std = p.nan
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname,'GeneralData.dat'):
                data = p.genfromtxt(filepath)
                if data[-1,4] != 0.0:
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
                        mean_turbul[j, 1] = data[k - frame_half_size : k \
                            + frame_half_size, 1].var()
                        k = k + 1
                    Grand_var = mean_turbul[:, 1].mean()
                    Grand_std = mean_turbul[:, 1].std()
                    Super_total_var = data[:, 1].var()
                    arg.append((Grand_var, Grand_std, turb_param,
                                rand_death_factor, Super_total_var))
                else:
                    break

GrandMeansData = []
os.path.walk(os.getcwd(), LoadEnvCond, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])

TheMeans = p.zeros((len(GrandMeansData), 4))
i = 0
for item in GrandMeansData:
    TheMeans[i,0] = item[0]
    TheMeans[i,1] = item[1]
    TheMeans[i,2] = item[2]
    TheMeans[i,3] = item[4]
    i += 1
xx = '# %(num)1.6f'%{"num":GrandMeansData[0][3]}
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)] # removing NaN's from the output
p.savetxt("VarianceOfEnvEtAl.dat", TheMeans, fmt='%1.6f')
doc = open('VarianceOfEnvEtAl.dat', 'a')
doc.write(xx)

print TheMeans

p.figure(9, figsize=(1000,600))
p.plot(TheMeans[:, 2], TheMeans[:, 0], 'ko--')
p.fill_between(TheMeans[:, 2], TheMeans[:, 0] + TheMeans[:, 1],
               TheMeans[:, 0] - TheMeans[:, 1], color=(0.75, 0.75, 0.75, 0.75))
p.plot(TheMeans[:, 2], TheMeans[:, 3], 'ko-')
p.xlabel('turbulence level ($ T $)', fontsize=LabelFontSize)
p.ylabel("variance of the env conditions",
         fontsize=LabelFontSize)
p.ylim(ymin=0)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.axis([0, 0.5, 0, 0.40001])
p.grid(True)
p.show()
