#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots averaged population age distribution in a form of histogram with
theoretical exponentiate distribution of cells age if it would depend only on
random death rate $ \delta $

Created, on, Mon, Oct, 11, 02:53:53, 2010

Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import re
import linecache as ln

LabelFontSize = 26
TickSize = 22
titleSize = 26

par_0 = re.split(" ", ln.getline('ModelParams.dat', 29))
death_rate = float(par_0[6])
age_data = p.genfromtxt("CellsAgeData.dat")
means = p.zeros((age_data.shape[1]))
STD = p.zeros((age_data.shape[1]))
age = p.arange(10, 10*age_data.shape[1]+10, 10)
summ = p.zeros((age_data.shape[0]))

age_data = age_data[p.ceil(age_data.shape[0]/2):, :]
for j in xrange(age_data.shape[0]):
#    summ[j] = age_data[j,:].sum()
    age_data[j, :] = age_data[j, :] / age_data[j, :].sum()


for i in xrange(age_data.shape[1]):
    means[i] = age_data[:, i].mean()
    STD[i] = age_data[:, i].std()

exp_distr = 10.0 * death_rate * p.exp(-death_rate * age)
# why the hell I need to multiply by 10 !? YES! 10 is the width of the
# histrogram bin!

maxAge = age.max() - 100
annot = "C"

p.figure(14, figsize=(12, 8))
p.bar(age, means, width=10.0, color='k', linewidth=0)
p.fill_between(age, means + STD, means - STD, color=(0.75, 0.75, 0.75, 0.75))
p.plot(age, exp_distr, 'w--', linewidth=3)
p.axis([0, maxAge, 0, 0.1])
p.xlabel('age of cells (steps)', fontsize=LabelFontSize)
p.ylabel('frequency of occurrence', fontsize=LabelFontSize)
p.title('$ \delta $ = %s' % death_rate, fontsize=titleSize)
p.text(1300, 0.09, annot, horizontalalignment='center',
       verticalalignment='center', fontsize=35)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)

print "Check sums"
print 'Means sum :', means.sum()
print 'Exp. sum  :', exp_distr.sum()

p.show()
