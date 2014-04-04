#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates mean surface of genotype's envelope +/- STD in number of
runs. Each run's results have to be stored in separate directory.
Produces a file with those results.
Change the value of 'important_line' variable to extract desirable parameter.

Created on Thu Oct  7 12:21:21 2010
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import os
import re
import sys
import linecache as ln
"""
data = genfromtxt("GeneralData.dat")
envelopes_sizes = genfromtxt("FrameSizeData.dat")
envelopes_maximums = genfromtxt("FrameMaxData.dat")
genome_size = genfromtxt("GenomeSizeData.dat")
age_data = genfromtxt("CellsAgeData.dat")
"""

try:
    important_line = int(sys.argv[1])
except:
    important_line = 8


def LoadEnvelopeMeanSize(arg, dirname, files):
    for file in files:
        Grand_mean = p.nan
        Grand_STD = p.nan
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GeneralData.dat'):
            data = p.genfromtxt(filepath)
            if data[-1, 4] != 0.0:
                data_chopped = data[1000:-1, :]
                Grand_mean = data_chopped[:, 12].mean()
                Grand_STD = p.sqrt((sum(data_chopped[:, 4]
                    * data_chopped[:, 13]**2) + sum((data_chopped[:, 12]
                    - Grand_mean)**2)) / sum(data_chopped[:, 4]))
                filepath_turb = os.path.join(dirname, 'ModelParams.dat')
                l = re.split(" ", ln.getline(filepath_turb, 6))
                turb_param = float(l[6])
            # change the line number in ln.getline below to extract the
            # parameter you're interested in
                l = re.split(" ", ln.getline(filepath_turb, important_line))
                interesting_param = float(l[6])
                arg.append((Grand_mean, Grand_STD, turb_param,
                            interesting_param))
            else:
                break


GrandMeansData = []
os.path.walk(os.getcwd(), LoadEnvelopeMeanSize, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])

#  ---- loading and checking data ---
TheMeans = p.zeros((len(GrandMeansData), 3))
i = 0
for item in GrandMeansData:
    TheMeans[i, 0] = item[0]
    TheMeans[i, 1] = item[1]
    TheMeans[i, 2] = item[2]
    i += 1
xx = '# %(num)1.6f' % {"num": GrandMeansData[0][3]}

print TheMeans
print xx
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)]  # removing NaN from the output

p.savetxt("Evelopes.dat", TheMeans, fmt='%1.6f')
doc = open('Evelopes.dat', 'a')
doc.write(xx)

p.figure(9, figsize=(10, 6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.plot(TheMeans[:, 2], TheMeans[:, 0], 'ko-')
#p.plot(TheMeans[:,2],TheMeans[:,0]-TheMeans[:,1], 'k--' )
#p.plot(TheMeans[:,2],TheMeans[:,0]+TheMeans[:,1], 'k--' )
#p.errorbar(TheMeans[:,2],TheMeans[:,0], yerr=TheMeans[:,1], fmt='o-', ecolor='k')
p.fill_between(TheMeans[:, 2], TheMeans[:, 0] + TheMeans[:, 1],
               TheMeans[:, 0] - TheMeans[:, 1], color=(0.75, 0.75, 0.75, 0.75))
p.axis([0, 0.5, 0, (TheMeans[:, 0] + TheMeans[:, 1]).max()])
p.xlabel('Turbulence level ($ T $)', fontsize=16)
p.ylabel('Genotype\'s envelope size', fontsize=16)
p.grid(True)
p.show()
