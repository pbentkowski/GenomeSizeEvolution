#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:29:34 2010

@author: piotr bentkowski
p.bentkowski@uea.ac.uk
"""
import pylab as p
import matplotlib.cm as cm

start = 7000
end = 10000
# start:end
# ---- importing and calculating data ---
data = p.genfromtxt("GeneralData.dat")
genome_size = p.genfromtxt("GenomeSizeData.dat")
age_data = p.genfromtxt("CellsAgeData.dat")
rel_genome_size = p.zeros((genome_size.shape[0],genome_size.shape[1]))
rel_age_data = p.zeros((age_data.shape[0],age_data.shape[1]))

# calculating stuff for genome size distribution plot
for i in range(genome_size.shape[0]):
        S = sum(genome_size[i, :])
        for j in range(genome_size.shape[1]):
            rel_genome_size[i,j] = genome_size[i, j] / S

# calculating stuff for age of cells distribution plot
for k in range(age_data.shape[0]):
        AgeSum = sum(age_data[k, :])
        for l in range(age_data.shape[1]):
            rel_age_data[k,l] = age_data[k,l] / AgeSum

frame_half_size = 300
mean_shannon = p.zeros((data.shape[0]-(2*frame_half_size),3))
mean_shannon[:,0] = data[frame_half_size:-frame_half_size][:,0].copy()
k = frame_half_size

for j in range(mean_shannon.shape[0]):
        mean_shannon[j,1] = data[k-frame_half_size:k+frame_half_size,9].mean()
        mean_shannon[j,2] = data[k-frame_half_size:k+frame_half_size,9].std()
        k = k + 1;

# ---- end of importing and calculating data  ---

#---- first plot-------
p.figure(1, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+10+10")

p.subplot(311)
p.plot(data[:, 0], data[:, 1],'k-')
p.axis([(data[start:end,0]).min(), (data[start:end,0]).max(), -1, 1])
p.ylabel('environmental conditions', fontsize=16)
p.grid(True)

p.subplot(312)
p.plot(data[:, 0],data[:, 4],'k-')
p.axis([(data[start:end,0]).min(), (data[start:end,0]).max(), 0, data[:,4].max()])
p.ylabel('number of cells', fontsize=16)
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 2],'k-', linewidth=1)
p.fill_between(data[:, 0],data[:, 2] + data[:, 3],data[:, 2] - data[:, 3], color=(0.75,0.75,0.75,0.75))
p.axis([(data[start:end,0]).min(), (data[start:end,0]).max(), 0, (data[:, 2] + data[:, 3]).max()])
p.xlabel('time (steps)', fontsize=16)
p.ylabel('mean number of genes', fontsize=16)
p.grid(True)
#---- end of first plot-------

#---- secont plot-------
p.figure(2, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+100+60")

p.subplot(311)
p.plot(data[:, 0], data[:, 5], 'k-')
p.axis([data[start:end, 0].min(),  data[start:end, 0].max(), 0,  data[start:end, 5].max()  ])
p.ylabel('resources inside cells', fontsize=16)
p.grid(True)

p.subplot(312)
p.plot(data[:, 0],data[:, 6], 'k-')
p.axis([data[start:end, 0].min(),  data[start:end, 0].max(), 0,  data[start:end, 6].max()])
p.ylabel('resources in the environment', fontsize=16)
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 5]+data[:, 6], 'k-')
p.axis([data[start:end, 0].min(),  data[start:end, 0].max(), (data[start:end, 5]+data[start:end, 6]).min(),  (data[start:end, 5]+data[start:end, 6]).max()])
p.xlabel('time (steps)', fontsize=16)
p.ylabel('diff. res.', fontsize=16)
p.grid(True)
#---- end of second plot-------

#---- third plot-------
x = data[:, 0]
y = p.r_[1:genome_size.shape[1]+1]
p.figure(3, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.contourf(x, y, rel_genome_size.transpose(), cmap=cm.Greys)
p.axis([data[start:end, 0].min(),  data[start:end, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)', fontsize=20)
p.ylabel('number of genes', fontsize=20)
p.xticks(size=16)
p.yticks(size=16)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence', fontsize=20)
#---- end of third plot-------


#---- fourth plot-------
p.figure(4, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+300+160")

p.subplot(311)
p.plot(data[:, 0], data[:, 11],'k-')
p.axis([data[start:end, 0].min(), data[start:end, 0].max(), 0.9*data[start:end, 11].min(),  1.1*data[start:end, 11].max()])
p.ylabel('number of dead', fontsize=16)
p.grid(True)

p.subplot(312)
p.plot(data[:, 0], data[:, 10],'k-')
p.axis([data[:, 0].min(),  data[start:end, 0].max(), 0.9*data[start:end, 10].min(),  1.1*data[start:end, 10].max()])
p.ylabel('number of born', fontsize=16)
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 9], 'k-')
p.plot(mean_shannon[:,0], mean_shannon[:,1],'k-', label='stepping mean',linewidth=2)
p.axis([data[start:end, 0].min(),  data[start:end, 0].max(), 0.9*data[start:end, 9].min(),  1.1*data[start:end, 9].max()])
p.xlabel('time (steps)', fontsize=16)
p.ylabel('Shannon index', fontsize=16)
p.grid(True)
#---- end of fourth plot-------

#---- fifth plot-------
yy = p.r_[1:age_data.shape[1]+1]
p.figure(5, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+400+210")
p.contourf(x, yy, p.log(rel_age_data.transpose()), cmap=cm.Greys)
#p.contourf(10.*x, yy, rel_age_data.transpose(), cmap=cm.Greys)
p.axis([data[start:end, 0].min(),  data[start:end, 0].max(),  yy.min(),  yy.max()])
p.xlabel('time (steps)', fontsize=16)
p.ylabel('Age of cells (steps)', fontsize=16)
p.grid(True)
cbb = p.colorbar(format=r"%.2f")
cbb.set_label('Logarithm of frequency of occurrence', fontsize=16)
#---- end of fifth plot-------

p.show()
print "done!"
