#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This code was design just to make a plot in the thesis. It the plot 4.5 titled:
Mean genome size and span of genome sizes grows with rising HGT probability.

Created on Fri Sep 28 21:14:03 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import pylab as p
import os
import re
import linecache as ln

AxisLabelFontSize = 17
AxisTickFontSize = 17
AnnotateFontSize = 20
interest_line = 53

#--- read data from files in directiories ---
def LoadGenomeMeanSize(arg, dirname, files):
       for file in files:
            Grand_mean = p.nan
            Grand_STD = p.nan
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname,'GeneralData.dat'):
                data = p.genfromtxt(filepath)
                if data[-1,4] != 0.0:
                    data_chopped = data[1000:-1,:]
                    Grand_mean = data_chopped[:,2].mean()
                    Grand_STD = p.sqrt((sum(data_chopped[:, 4] \
                    * data_chopped[:, 3]**2) + sum((data_chopped[:, 2] \
                    - Grand_mean)**2)) / sum(data_chopped[:, 4]))
                    filepath_hgt = os.path.join(dirname,'ModelParams.dat')
                    l = re.split(" ", ln.getline(filepath_hgt, interest_line))
                    hgt_param = float(l[6])
                    l_1 = re.split(" ", ln.getline(filepath_hgt, 6))
                    T_value = float(l_1[6])
                    arg.append((Grand_mean, Grand_STD, hgt_param, T_value))
                else:
                    break

GrandMeansData = []
os.path.walk(os.getcwd(), LoadGenomeMeanSize, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])

#  ---- loading and cheking data ---
TheMeans = p.zeros((len(GrandMeansData), 4))
i = 0
for item in GrandMeansData:
    TheMeans[i, 0] = item[0]
    TheMeans[i, 1] = item[1]
    TheMeans[i, 2] = item[2]
    TheMeans[i, 3] = item[3]
    i += 1

print TheMeans
TheMeans = TheMeans[~p.isnan(TheMeans).any(1)] # removing NaN's from the output

p.savetxt("GenomeSizes.dat", TheMeans, fmt='%1.6f')

XX = 0.00125
YY = 36
p.figure(10, figsize=(10,6))
#p.subplot(221)
p.plot(TheMeans[:,2],TheMeans[:,0], 'ko-' )
p.fill_between(TheMeans[:,2], TheMeans[:,0] + TheMeans[:,1],
               TheMeans[:,0] - TheMeans[:,1], color=(0.75,0.75,0.75,0.75))
p.annotate('$T = %1.2f$'%(TheMeans[0,3],), xy = (XX, YY),  xycoords='data',
           fontsize=AnnotateFontSize)
p.ylim( 0, 40)
p.xlabel('HGT probability on cell level ($ h_{c} $)', fontsize=AxisLabelFontSize)
p.ylabel('number of genes', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
p.show()