#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 12:20:21 2010

@author: Piotr Bentkowski
"""

import os
import re
import pylab as p
import numpy as np
#from StringIO import StringIO
import linecache as ln
#import string as str
#import matplotlib.cm as cm


#--- read data from files in directiories ---
def LoadMyData(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'GeneralData.dat'):
                data = np.genfromtxt(filepath)
                Shannon = data[:, 9]
            if filepath == os.path.join(dirname, 'ModelParams.dat'):
                l = re.split(" ", ln.getline(filepath, 6))
                turb_param = float(l[6])
                arg.append((Shannon, turb_param))

AllData = []
os.path.walk(os.getcwd(), LoadMyData, AllData)
AllData = sorted(AllData, key=lambda data_sort: data_sort[1])

#  ---- loading and cheking data ---
SHANNON = p.zeros((len(AllData[0][0]), len(AllData)))
tags = p.zeros((len(AllData), 1))
i = 0
temp_item_size = p.size(p.array(AllData[0][0]))
for item in AllData:
        item_size = p.size(item[0])
        if temp_item_size != item_size:
            print 'ERROR IN DATA: At least one of the data series has',
            print 'a different length! Procedure aborted. Check your data and',
            print 'try again or don\'t ever come back...'
            exit()
        else:
            SHANNON[:, i] = item[0]
            tags[i] = item[1]
            i = i + 1
            temp_item_size = item_size
print 'This job has', SHANNON.shape[1], 'data series with', SHANNON.shape[0], \
      'data point each.'

frame_min = 2  # or other even number
frame_max = 200
frame_step = 2  # or other even number
if frame_min % 2 == 1:   # frame size has to be an even number
    frame_min = frame_min + 1
if frame_max % 2 == 1:   # frame size has to be an even number
    frame_max = frame_max - 1

All_the_frames = p.arange(frame_min, frame_max, frame_step)
Mean_Var = p.zeros((All_the_frames.shape[0], SHANNON.shape[1]))

# --- caunting the variances ---
for i in range(SHANNON.shape[1]):
        print 'Processing data series No.', i+1, 'out of', SHANNON.shape[1]
        l = 0
        for item in All_the_frames:
            diff_of_means_size = SHANNON.shape[0] - item
            diff_of_means = p.zeros((diff_of_means_size, 1))
            half = item/2
            k = half
            for j in range(diff_of_means.shape[0]):
                diff_of_means[j] = SHANNON[k-half:k+half, i].var()
                k = k + 1
            Mean_Var[l, i] = diff_of_means.sum() / float(diff_of_means_size)
            l = l + 1
        i = i + 1
print 'I\'m done with processing data series. Writing to files.'

np.savetxt("SHN_variance.dat", Mean_Var)
np.savetxt("SHN_TAGS_variance.dat", tags)
np.savetxt("SHN_All_the_frames.dat", All_the_frames)

#p.figure(1)
#for i in range(SHANNON.shape[1]):
#        p.plot(All_the_frames, Mean_Var[:,i], '-', label='%4.2f'%(tags[i],),
#                linewidth=2)
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
#p.legend(fancybox=True)
##p.axis([0.0,  frame_max,  0.0,  1.1*Mean_Var.max()])
#p.xlabel('frame size')
#p.ylabel('mean variance of the Shannon index')
#p.grid(True)
#p.show()

print 'Over. Have a a nice plot! :-)'
