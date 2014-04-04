#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Mon Mar 17 12:49:30 2014
Author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""
import pylab as p
import matplotlib.cm as cm

# ---- importing and calculating data ---
data = p.genfromtxt("GeneralData.dat")
genome_size = p.genfromtxt("GenomeSizeData.dat")
age_data = p.genfromtxt("CellsAgeData.dat")
rel_genome_size = p.zeros((genome_size.shape[0], genome_size.shape[1]))
rel_age_data = p.zeros((age_data.shape[0], age_data.shape[1]))

# calculating stuff for genome size distribution plot
for i in range(genome_size.shape[0]):
    S = sum(genome_size[i, :])
    for j in range(genome_size.shape[1]):
        rel_genome_size[i, j] = genome_size[i, j] / S

# calculating stuff for age of cells distribution plot
for k in range(age_data.shape[0]):
    AgeSum = sum(age_data[k, :])
    for l in range(age_data.shape[1]):
        rel_age_data[k, l] = age_data[k, l] / AgeSum

frame_half_size = 300
mean_shannon = p.zeros((data.shape[0]-(2*frame_half_size), 3))
mean_shannon[:, 0] = data[frame_half_size:-frame_half_size][:, 0].copy()
k = frame_half_size

for j in range(mean_shannon.shape[0]):
    mean_shannon[j, 1] = data[k-frame_half_size:k+frame_half_size, 9].mean()
    mean_shannon[j, 2] = data[k-frame_half_size:k+frame_half_size, 9].std()
    k = k + 1

# ---- end of importing and calculating data  ---

#--- Fonts sizes ---
X_ticks = 23
Y_ticks = 23#!/usr/bin/env python
import pylab as p
import matplotlib.cm as cm

# ---- importing and calculating data ---
data = p.genfromtxt("GeneralData.dat")
genome_size = p.genfromtxt("GenomeSizeData.dat")
age_data = p.genfromtxt("CellsAgeData.dat")
rel_genome_size = p.zeros((genome_size.shape[0], genome_size.shape[1]))
rel_age_data = p.zeros((age_data.shape[0], age_data.shape[1]))

# calculating stuff for genome size distribution plot
for i in range(genome_size.shape[0]):
    S = sum(genome_size[i, :])
    for j in range(genome_size.shape[1]):
        rel_genome_size[i, j] = genome_size[i, j] / S

# calculating stuff for age of cells distribution plot
for k in range(age_data.shape[0]):
    AgeSum = sum(age_data[k, :])
    for l in range(age_data.shape[1]):
        rel_age_data[k, l] = age_data[k, l] / AgeSum

frame_half_size = 300
mean_shannon = p.zeros((data.shape[0]-(2*frame_half_size), 3))
mean_shannon[:, 0] = data[frame_half_size:-frame_half_size][:, 0].copy()
k = frame_half_size

for j in range(mean_shannon.shape[0]):
    mean_shannon[j, 1] = data[k-frame_half_size:k+frame_half_size, 9].mean()
    mean_shannon[j, 2] = data[k-frame_half_size:k+frame_half_size, 9].std()
    k = k + 1

# ---- end of importing and calculating data  ---

#--- Fonts sizes ---
X_ticks = 23
Y_ticks = 23
X_Label = 25
Y_Label = 25
Figs = (20, 12)
#---- first plot-------
p.figure(1, figsize=Figs)

p.subplot(311)
p.plot(data[:, 0], data[:, 1], 'k-')
p.axis([0, data[:, 0].max(), -1, 1])
p.ylabel('environmental\n conditions', fontsize=X_Label)
p.xticks(size=X_ticks)
p.yticks(size=Y_ticks)
p.grid(True)

p.subplot(312)
p.plot(data[:, 0], data[:, 4], 'k-')
p.axis([0, data[:, 0].max(), 0, data[:, 4].max()])
p.ylabel('number of cells', fontsize=X_Label)
p.xticks(size=X_ticks)
p.yticks(size=Y_ticks)
p.grid(True)

p.subplot(313)
x = data[:, 0]
y = p.r_[1:genome_size.shape[1]+1]
p.contourf(x, y, rel_genome_size.transpose(), cmap=cm.binary)
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)', fontsize=X_Label)
p.ylabel('number of genes', fontsize=Y_Label)
p.xticks(size=X_ticks)
p.yticks(size=Y_ticks)
p.grid(True)
cb = p.colorbar(format=r"%.2f", orientation='horizontal', aspect=40)
cb.ax.tick_params(labelsize=X_ticks)
cb.set_label('frequency of occurrence', fontsize=X_Label)
p.show()
