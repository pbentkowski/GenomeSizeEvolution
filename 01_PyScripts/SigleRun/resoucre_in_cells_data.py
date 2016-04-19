#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Takes files RecourceInCells.dat and RecourceUptakenCells.dat and show histograms
of average uptake values and distribution of resource allocation among the cell
for the stable second half of the model iteration.

Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import re
import linecache as ln

LabelsFontSize = 22
AxisTickFontSize = 16

F2_LabelsFontSize = 22
F2_AxisTickFontSize = 20

par_0 = re.split(" ", ln.getline('ModelParams.dat', 6))
tubulence = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 23))
max_intake_value = float(par_0[6])
par_1 = re.split(" ", ln.getline('ModelParams.dat', 25))
death_value = float(par_1[6])
par_2 = re.split(" ", ln.getline('ModelParams.dat', 26))
reprod_value = float(par_2[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 29))
death_rate = float(par_0[6])
l_3 = re.split(" ", ln.getline('ModelParams.dat', 37))
l_2 = re.split(" ", ln.getline('ModelParams.dat', 38))
l_1 = re.split(" ", ln.getline('ModelParams.dat', 39))
hist_beginning = float(l_1[6])
hist_bin_size = float(l_2[6])
hist_end = hist_beginning + float(l_3[6]) * hist_bin_size

inside_res_data = p.genfromtxt("RecourceInCells.dat")
half_size = p.floor(inside_res_data.shape[0] / 2)
inside_res_data = inside_res_data[half_size:-1, :]
means_in = p.zeros((inside_res_data.shape[1]))
STD_in = p.zeros((inside_res_data.shape[1]))
ranges = p.arange(hist_beginning, hist_end, hist_bin_size)
summ = p.zeros((inside_res_data.shape[0]))

uptake_res_data = p.genfromtxt("RecourceUptakenCells.dat")
uptake_res_data = uptake_res_data[half_size:-1, :]

means_uptake = p.zeros((uptake_res_data.shape[1]))
STD_uptake = p.zeros((uptake_res_data.shape[1]))
ranges_uptake = p.arange(0.0, uptake_res_data.shape[1], 1.0)
summ = p.zeros((uptake_res_data.shape[0]))

for j in xrange(inside_res_data.shape[0]):
    inside_res_data[j,:] = inside_res_data[j,:] / inside_res_data[j,:].sum()

for i in xrange(inside_res_data.shape[1]):
    means_in[i] = inside_res_data[:,i].mean()
    STD_in[i] = inside_res_data[:,i].std()

for k in xrange(uptake_res_data.shape[0]):
    uptake_res_data[k,:] = uptake_res_data[k,:] / uptake_res_data[k,:].sum()

for l in xrange(uptake_res_data.shape[1]):
    means_uptake[l] = uptake_res_data[:,l].mean()
    STD_uptake[l] = uptake_res_data[:,l].std()

hungry_ones = means_uptake[0]
STD_of_hungry = STD_uptake[0]
means_uptake[0] = 0.0
STD_uptake[0] = 0.0

X_max = p.ceil((means_uptake + STD_uptake).max() * 100.0) / 100.0
Y_max = p.ceil((means_in + STD_in).max() * 100.0) / 100.0
p.figure(1, figsize=(20, 14), dpi=80)
p.subplot(221)
p.bar(ranges, means_in + STD_in, width = hist_bin_size,
         color=(0.75, 0.75, 0.75, 0.75), linewidth=0)
p.bar(ranges, means_in, width = hist_bin_size, color='k', linewidth=0)
p.vlines(death_value, 0.0, Y_max, 'k', linestyles='dashed')
p.vlines(reprod_value, 0.0, Y_max, 'k', linestyles='dashed')
#p.fill_between(ranges, means_in + STD_in, means_in - STD_in,
#                color=(0.75, 0.75, 0.75, 0.75))
p.axis([ranges.min(), ranges.max(), 0.0, Y_max])
p.xlabel('resources possessed by the cells', fontsize=LabelsFontSize)
p.ylabel('frequency of occurrence', fontsize=LabelsFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
p.title(r'$\delta = %(death)1.3f$   $T= %(turb)1.2f$'%{"death":death_rate,
        "turb":tubulence}, fontsize=LabelsFontSize)

#p.figure(2, figsize=(1000, 600))
p.subplot(223)
#p.get_current_fig_manager().window.wm_geometry("1000x600+10+10")
p.bar(ranges_uptake, means_uptake + STD_uptake, width = 1.0,
         color=(0.75, 0.75, 0.75, 0.75), linewidth=0)
p.bar(ranges_uptake, means_uptake, width = 1.0, color='k', linewidth=0)
p.vlines(max_intake_value, 0.0, X_max, 'k', linestyles='dashed')
#p.fill_between(ranges_uptake, means_uptake + STD_uptake,
#                means_uptake - STD_uptake, color=(0.75, 0.75, 0.75, 0.75))
p.axis([ranges_uptake.min(), ranges_uptake.max(), 0.0, X_max])
p.xlabel('resources up-taken by the cells', fontsize=LabelsFontSize)
p.ylabel('frequency of occurrence', fontsize=LabelsFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
#p.title(r'$\delta = %(death)1.3f$   $T= %(turb)1.2f$'%{"death":death_rate,
#        "turb":tubulence}, fontsize=LabelsFontSize)

#p.savefig('Fig1.png',bbox_inches='tight', pad_inches=0.1)

print ' Total fraction of collected resources (just cheking):', means_in.sum()
print 'Fraction of hungry cells:', hungry_ones, '+/-', STD_of_hungry


p.figure(2, figsize=(10, 6))
p.bar(ranges_uptake, means_uptake + STD_uptake, width = 1.0,
         color=(0.75, 0.75, 0.75, 0.75), linewidth=0)
p.bar(ranges_uptake, means_uptake, width = 1.0, color='k', linewidth=0)
p.vlines(max_intake_value, 0.0, X_max, 'k', linestyles='dashed')
#p.fill_between(ranges_uptake, means_uptake + STD_uptake,
#                means_uptake - STD_uptake, color=(0.75, 0.75, 0.75, 0.75))
p.axis([ranges_uptake.min(), ranges_uptake.max(), 0.0, X_max])
p.xlabel('resources up-taken by the cells', fontsize=F2_LabelsFontSize)
p.ylabel('frequency of occurrence', fontsize=F2_LabelsFontSize)
p.xticks(size=F2_AxisTickFontSize)
p.yticks(size=F2_AxisTickFontSize)
p.grid(True)
p.title(r'$T= %(turb)1.2f$'% {"turb":tubulence}, fontsize=F2_LabelsFontSize)

p.show()