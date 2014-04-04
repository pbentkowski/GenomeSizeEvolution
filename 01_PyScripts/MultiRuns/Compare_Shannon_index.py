#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the anarage Shannon index for model's runs located in different
directories. Change the value of 'important_line' variable to extract desirable
parameter.

Created on Mon May 30 21:34:51 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import sys
import pylab as p
import os
import re
import linecache as ln

important_line = int(sys.argv[1])


#--- read data from files in directiories ---
def LoadShannonIndex(arg, dirname, files):
    for file in files:
        Grand_mean = p.nan
        Grand_STD = p.nan
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GeneralData.dat'):
            data = p.genfromtxt(filepath)
            if data[-1, 4] != 0.0:
                half_data = p.ceil(data.shape[0] / 2.0)
                data_chopped = data[half_data:-1, :]
                Grand_mean = data_chopped[:, 9].mean()
                Grand_STD = data_chopped[:, 9].std()
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
os.path.walk(os.getcwd(), LoadShannonIndex, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])

#  ---- loading and checking data ---
TheMeans = p.zeros((len(GrandMeansData), 3))
i = 0
for item in GrandMeansData:
    TheMeans[i, 0] = item[0]
    TheMeans[i, 1] = item[1]
    TheMeans[i, 2] = item[2]
    i += 1
xx = '# %(num)1.4f' % {"num": GrandMeansData[0][3]}

print TheMeans
print xx
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)]  # removing NaN from the output

p.savetxt("Shannons.dat", TheMeans, fmt='%1.5f')
doc = open('Shannons.dat', 'a')
doc.write(xx)

p.figure(19, figsize=(12, 6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.plot(TheMeans[:, 2], TheMeans[:, 0], 'ko-')
#p.plot(TheMeans[:,2],TheMeans[:,0]-TheMeans[:,1], 'k--' )
#p.plot(TheMeans[:,2],TheMeans[:,0]+TheMeans[:,1], 'k--' )
#p.errorbar(TheMeans[:,2],TheMeans[:,0], yerr=TheMeans[:,1], fmt='o-', ecolor='k')
p.fill_between(TheMeans[:, 2], TheMeans[:, 0] + TheMeans[:, 1],
               TheMeans[:, 0] - TheMeans[:, 1], color=(0.75, 0.75, 0.75, 0.75))
p.axis([0, 0.5, 0, (TheMeans[:, 0] + TheMeans[:, 1]).max()])
p.xlabel('Turbulence level ($ T $)', fontsize=16)
p.ylabel('Population\'s mean Shannon index', fontsize=16)
p.grid(True)
p.show()
