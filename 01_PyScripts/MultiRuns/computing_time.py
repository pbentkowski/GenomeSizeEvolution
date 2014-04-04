#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the average computing time using a bunch of model's data files stored
in different directories.

Created on Wed Sep 29 17:15:45 2010
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import os, re
import pylab as p
import linecache as ln
# the one below is mine
import file_len as fl

#--- read data from files in directiories ---
def LoadMyData(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname,'GeneralData.dat'):
                line_num = fl.file_len(filepath)
                l_1 = re.split(" ", ln.getline(filepath, line_num))
                time_of_comp = float(l_1[5])
                l_2 = re.split(" ", ln.getline(filepath, line_num - 1))
                how_many_cells = float(l_2[4])
                if how_many_cells != 0:
                   did_make_it = 1
                else:
                    did_make_it = 0
                arg.append((time_of_comp, did_make_it))

AllData = []
os.path.walk(os.getcwd(), LoadMyData, AllData)

ENV = p.zeros(( len(AllData), 2) )

i = 0
for item in AllData:
        ENV[i, 0] = item[0]
        ENV[i, 1] = item[1]
        i = i + 1

MeanTime = ENV[:, 0].mean()
STDtime = ENV[:, 0].std()

hours = p.floor(MeanTime/3600.0)
minutes = p.floor((MeanTime - hours * 3600.0) / 60.0)
seconds = p.floor((MeanTime - hours * 3600.0 - minutes * 60.0))

print "There were", len(AllData), "runs and", ENV[:, 1].sum(), "populations "\
"made it to the end of simulation."
print "The mean computing time was:"
print MeanTime, "+/-", STDtime, "seconds"
#print "That is", hours, "hours", minutes, "minutes and",  seconds, "seconds"
print 'That is %(hours)2.0f:%(minutes)2.0f:%(seconds)2.0f hours on avarege'% \
       {"hours": hours, "minutes": minutes, "seconds":seconds}