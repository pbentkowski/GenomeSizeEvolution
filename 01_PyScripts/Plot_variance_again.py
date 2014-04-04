#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 16:33:22 2010

@author: Piotr Bentkowski
"""
import pylab as p
import matplotlib.cm as cm

Mean_Var = p.genfromtxt("ENV_variance.dat")
tags = p.genfromtxt("ENV_TAGS_variance.dat")
All_the_frames = p.genfromtxt("ENV_All_the_frames.dat")


p.figure(3, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.contourf(All_the_frames, tags, Mean_Var.transpose(), cmap=cm.gist_yarg)
p.axis([All_the_frames.min(), All_the_frames.max(), 0, tags.max()])
p.xlabel('cell age (steps)', fontsize=16)
p.ylabel('turbulance level $ T $', fontsize=16)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('mean variance', fontsize=16)
p.show()