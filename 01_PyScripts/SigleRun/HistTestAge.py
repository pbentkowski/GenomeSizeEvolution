#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created, on, Mon, Oct, 11, 02:53:53, 2010

Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import re
import linecache as ln

LabelFontSize=25
TickSize=21

par_0 = re.split(" ", ln.getline('ModelParams.dat', 29))
death_rate = float(par_0[6])
age_data = p.genfromtxt("CellsAgeData.dat")
means = p.zeros((age_data.shape[1]))
STD = p.zeros((age_data.shape[1]))
age = p.arange(10, 10*age_data.shape[1]+10, 10 )
summ = p.zeros((age_data.shape[0]))

age_data = age_data[p.ceil(age_data.shape[0]/2):,:]
for j in xrange(age_data.shape[0]):
#    summ[j] = age_data[j,:].sum()
    age_data[j,:] = age_data[j,:] / age_data[j,:].sum()
    

for i in xrange(age_data.shape[1]):
    means[i] = age_data[:,i].mean()
    STD[i] = age_data[:,i].std()

exp_distr = 10.0 * death_rate * p.exp(-death_rate * age)
# why the hell I need to multiply by 10 !? YES! 10 is the width of the
# histrogram bin!

p.figure(14, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+10+10")
p.bar(age, means, width=10.0, color='k', linewidth=0)
p.fill_between(age, means + STD, means - STD, color=(0.75,0.75,0.75,0.75))
p.plot(age, exp_distr, 'w--', linewidth=3)
p.axis([0, age.max(), 0, 0.15])
#p.axis([0, age.max(), 0, (means + STD).max()])
p.xlabel('Age of cells (steps)', fontsize=LabelFontSize)
p.ylabel('Frequency of occurrence', fontsize=LabelFontSize)
p.title(death_rate, fontsize=27)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)

print '  ', means.sum()
print '  ', exp_distr.sum()

p.show()