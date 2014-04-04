#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 15:26:21 2010

@author: piotr bentkowski
p.bentkowski@uea.ac.uk
"""
import pylab as p
import matplotlib.cm as cm
import re
import linecache as ln

LabelsFontSize = 17
AxisTickFontSize = 16
data = p.genfromtxt("GeneralData.dat")
uptake_res_data = p.genfromtxt("RecourceUptakenCells.dat")
inside_res_data = p.genfromtxt("RecourceInCells.dat")
rel_uptake_res_data = p.zeros((uptake_res_data.shape[0], uptake_res_data.shape[1]))
rel_inside_res_data = p.zeros((inside_res_data.shape[0], inside_res_data.shape[1]))

par_0 = re.split(" ", ln.getline('ModelParams.dat', 23))
max_intake_value = float(par_0[6])
par_1 = re.split(" ", ln.getline('ModelParams.dat', 25))
death_value = float(par_1[6])
par_2 = re.split(" ", ln.getline('ModelParams.dat', 26))
reprod_value = float(par_2[6])
l_3 = re.split(" ", ln.getline('ModelParams.dat', 37))
l_2 = re.split(" ", ln.getline('ModelParams.dat', 38))
l_1 = re.split(" ", ln.getline('ModelParams.dat', 39))
hist_beginning = float(l_1[6])
hist_bin_size = float(l_2[6])
hist_end = hist_beginning + float(l_3[6]) * hist_bin_size

# Most of the cells did not got any resource, so we exclude them from the analysis.
uptake_res_data[:, 0] = 0.0

# calculating stuff for resource uptake plot
for i in range(uptake_res_data.shape[0]):
        S = sum(uptake_res_data[i, :])
        for j in range(uptake_res_data.shape[1]):
            rel_uptake_res_data[i,j] = uptake_res_data[i, j] / S

# calculating stuff for resources in the cells plot
for k in range(inside_res_data.shape[0]):
        AgeSum = sum(inside_res_data[k, :])
        for l in range(inside_res_data.shape[1]):
            rel_inside_res_data[k, l] = inside_res_data[k,l] / AgeSum

x_length = len(data[:, 0])
#---- first plot-------
x = data[:, 0]
y = p.r_[1:uptake_res_data.shape[1]+1]
y_max_intake = p.ones(x_length)*max_intake_value
p.figure(1, figsize=(10, 6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.contourf(x, y, rel_uptake_res_data.transpose(), cmap=cm.jet)
p.plot(x, y_max_intake, 'y--')
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)', fontsize = LabelsFontSize)
p.ylabel('Resource taken by the cells', fontsize = LabelsFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence', fontsize = LabelsFontSize)
#---- end of first plot-------

#---- second plot-------
x = data[:, 0]
y = p.arange(hist_beginning,hist_end,hist_bin_size)
y_death = p.ones(x_length)*death_value
y_repr = p.ones(x_length)*reprod_value
p.figure(2, figsize=(10,6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.contourf(x, y, rel_inside_res_data.transpose(), cmap=cm.jet)
p.plot(x, y_death, 'r--')
p.plot(x, y_repr, 'g--')
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)', fontsize = LabelsFontSize)
p.ylabel('Resource inside the cells', fontsize = LabelsFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence', fontsize = LabelsFontSize)
#---- end of secont plot-------

p.show()
