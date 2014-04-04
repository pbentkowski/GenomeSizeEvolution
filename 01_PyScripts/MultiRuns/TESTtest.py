#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This scripts calculates mean number of genes +/- STD in the number of runs. 
Each run's results have to be stored in separate directory. 
Produces a file with those results.
"""

import pylab as p
import os, re
import linecache as ln

def LoadGenomeMeanSize(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname,'GeneralData.dat'):
                data = p.genfromtxt(filepath)
                if data[-1,4] != 0.0: # checking if data set is OK 
                    data_chopped = data[1000:-1,:]  # removing some of data
                    Grand_mean = data_chopped[:,2].mean()
                    Grand_STD = p.sqrt((sum(data_chopped[:,4] * data_chopped[:,3]**2)
                    + sum((data_chopped[:,2]-Grand_mean)**2))/sum(data_chopped[:,4]))

                    if filepath == os.path.join(dirname,'ModelParams.dat'):
                        l = re.split(" ", ln.getline(filepath, 6))
                        turb_param = float(l[2])                
                        arg.append((Grand_mean, Grand_STD, turb_param))
                else:
                    break

   
             
GrandMeansData = []
os.path.walk(os.getcwd(), LoadGenomeMeanSize, GrandMeansData)
GrandMeansData = sorted(GrandMeansData, key=lambda data_sort: data_sort[2])

TheMeans = p.zeros((len(GrandMeansData), 3 ))
i = 0
for item in GrandMeansData:
    TheMeans[i,0] = item[0]
    TheMeans[i,1] = item[1]
    TheMeans[i,2] = item[2]
    i = i + 1

print TheMeans # just checking...
# later do some computation on TheMeans in NumPy
