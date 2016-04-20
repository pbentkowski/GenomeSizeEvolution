#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 19:05:54 2010

@author: Piotr Bentkowski
"""

import os
import re
import pylab as p
import numpy as np
import linecache as ln


#--- read data from files in directiories ---
def LoadMyData(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'GeneralData.dat'):
                data = np.genfromtxt(filepath)
                Env = data[:, 1]
            if filepath == os.path.join(dirname, 'ModelParams.dat'):
                l = re.split(" ", ln.getline(filepath, 6))
                turb_param = float(l[2])
                arg.append((Env, turb_param))

AllData = []
os.path.walk(os.getcwd(), LoadMyData, AllData)
AllData = sorted(AllData, key=lambda data_sort: data_sort[1])

#  ---- loading and cheking data ---
ENV = p.zeros((len(AllData[0][0]), len(AllData)))
tags = p.zeros((len(AllData), 1))
i = 0
temp_item_size = p.size(p.array(AllData[0][0]))
for item in AllData:
        item_size = p.size(item[0])
        if temp_item_size != item_size:
            print 'ERROR IN DATA: At least one of the data series has a ',
            print 'different length! Procedure aborted. Check your data and',
            print 'try again or don\'t ever come back...'
            exit()
        else:
            ENV[:, i] = item[0]
            tags[i] = item[1]
            i = i + 1
            temp_item_size = item_size
print 'This job has', ENV.shape[1], 'data series with', ENV.shape[0], \
      'data point each.'

frame_min = 2  # or other even number
frame_max = 200
frame_step = 2  # or other even number
if frame_min % 2 == 1:   # frame size has to be an even number
    frame_min = frame_min + 1
if frame_max % 2 == 1:   # frame size has to be an even number
    frame_max = frame_max - 1

All_the_frames = p.arange(frame_min, frame_max, frame_step)
Mean_Var = p.zeros((All_the_frames.shape[0], ENV.shape[1]))

# --- caunting the variances ---
for i in xrange(ENV.shape[1]):
        print 'Processing data series No.', i+1, 'out of', ENV.shape[1]
        l = 0
        for item in All_the_frames:
            diff_of_means_size = ENV.shape[0] - item
            diff_of_means = p.zeros((diff_of_means_size, 1))
            half = item/2
            k = half
            for j in xrange(diff_of_means.shape[0]):
                diff_of_means[j] = ENV[k-half:k+half, i].var()
                k = k + 1
            Mean_Var[l, i] = diff_of_means.sum() / float(diff_of_means_size)
            l = l + 1
        i = i + 1
print 'I\'m done with processing data series. Writing to files.'

np.savetxt("ENV_variance.dat", Mean_Var)
np.savetxt("ENV_TAGS_variance.dat", tags)
np.savetxt("ENV_All_the_frames.dat", All_the_frames)

print 'Over. Have a a nice plot! :-)'
