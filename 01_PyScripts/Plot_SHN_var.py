#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 17:16:03 2010

@author: Piotr Bentkowski
"""
import matplotlib.cm as cm
import pylab as p
import numpy as np

Mean_Var = np.genfromtxt("SHN_variance.dat")
tags = np.genfromtxt("SHN_TAGS_variance.dat")
All_the_frames = np.genfromtxt("SHN_All_the_frames.dat")

p.figure(1, figsize=(1000,600))
for i in range(Mean_Var.shape[1]):
        ax = p.plot(All_the_frames, Mean_Var[:,i], '-', label='%4.2f'%(tags[i],), linewidth=2)
        temp_XX = int(round(All_the_frames.max()-(float(i)*(All_the_frames.max()/20.))))
        temp_YY = np.where(All_the_frames >= int(round(All_the_frames.max()-(float(i+1)*(All_the_frames.max()/20.)))))
        YY = Mean_Var[temp_YY[0][0],i]
        XX = All_the_frames[temp_YY[0][0]]
        ax = p.annotate('%4.2f'%(tags[i],), xy = (XX, YY),  xycoords='data')
##p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
#p.legend(fancybox=True)
#p.axis([0.0,  frame_max,  0.0,  1.1*Mean_Var.max()])
p.xlabel('frame size')
p.ylabel('mean variance of the Shannon index')
p.grid(True)
p.show()