#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This scripts calculates mean number of genes +/- STD in the number of runs.
Each run's results have to be stored in separate directory. Produces a file
with those results.

Created on Wed Oct  6 16:39:28 2010
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import os
import re
import linecache as ln
"""
data = genfromtxt("GeneralData.dat")
envelopes_sizes = genfromtxt("FrameSizeData.dat")
envelopes_maximums = genfromtxt("FrameMaxData.dat")
genome_size = genfromtxt("GenomeSizeData.dat")
age_data = genfromtxt("CellsAgeData.dat")
"""


def LoadGenomeMeanSize(arg, dirname, files):
    for file in files:
        Grand_mean = p.nan
        Grand_STD = p.nan
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GeneralData.dat'):
            data = p.genfromtxt(filepath)
            if data[-1, 4] != 0.0:
                data_chopped = data[1000:-1, :]
                Grand_mean = data_chopped[:, 2].mean()
                Grand_STD = p.sqrt((sum(data_chopped[:, 4]
                            * data_chopped[:, 3]**2) + sum((data_chopped[:, 2]
                            - Grand_mean)**2)) / sum(data_chopped[:, 4]))
                filepath_turb = os.path.join(dirname, 'ModelParams.dat')
                l = re.split(" ", ln.getline(filepath_turb, 6))
                turb_param = float(l[6])
                arg.append((Grand_mean, Grand_STD, turb_param))
            else:
                break

GrandMeansData = []
os.path.walk(os.getcwd(), LoadGenomeMeanSize, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])

#  ---- loading and cheking data ---
TheMeans = p.zeros((len(GrandMeansData), 3))
i = 0
for item in GrandMeansData:
    TheMeans[i, 0] = item[0]
    TheMeans[i, 1] = item[1]
    TheMeans[i, 2] = item[2]
    i += 1

print TheMeans
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)]  # removing NaN from the output

p.savetxt("GenomeSizes.dat", TheMeans, fmt='%1.6f')
p.figure(10, figsize=(10, 6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.plot(TheMeans[:, 2], TheMeans[:, 0], 'ko-')
#p.plot(TheMeans[:,2],TheMeans[:,0]-TheMeans[:,1], 'k--' )
#p.plot(TheMeans[:,2],TheMeans[:,0]+TheMeans[:,1], 'k--' )
#p.errorbar(TheMeans[:,2],TheMeans[:,0], yerr=TheMeans[:,1], fmt='o-', ecolor='k')
p.fill_between(TheMeans[:, 2], TheMeans[:, 0] + TheMeans[:, 1],
               TheMeans[:, 0] - TheMeans[:, 1], color=(0.75, 0.75, 0.75, 0.75))
#p.axis([0, 0.5, 2, 20])
p.xlabel('turbulence level ($ T $)', fontsize=16)
p.ylabel('number of genes', fontsize=16)
p.grid(True)
p.show()
