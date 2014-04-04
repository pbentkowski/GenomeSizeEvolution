# -*- coding: utf-8 -*-
"""
Plots general statistics of a run: env conditions, population size and average
gene number.

Created on Wed Jan 11 23:03:01 2012

Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import re
import linecache as ln

#--- Fonts sizes ---
X_ticks = 20
Y_ticks = 20
X_Label = 22
Y_Label = 22
TextFontSize = 24

data = p.genfromtxt("GeneralData.dat")
par_0 = re.split(" ", ln.getline('ModelParams.dat', 6))
tubulence = float(par_0[6])

p.figure(1, figsize=(18, 14))

p.subplot(311)
p.plot(data[:, 0], data[:, 1], 'k-')
p.axis([0, data[:, 0].max(), -1, 1])
p.ylabel('environmental\n conditions', fontsize=X_Label)
p.xticks(size=X_ticks)
p.yticks(size=Y_ticks)
p.grid(True)

p.subplot(312)
locality = data[-1, 0] - int(data[-1, 0] * 0.1)
tick_marks_01 = p.arange(0., 4000.0001, 1000.)
p.plot(data[:, 0], data[:, 4], 'k-')
p.axis([0, data[:, 0].max(), 0, data[:, 4].max()])
p.ylabel('number of cells', fontsize=X_Label)
p.xticks(size=X_ticks)
p.yticks(tick_marks_01, size=Y_ticks)
p.ylim(ymax=4000)
p.text(locality, 500.0, '$T = %(turb)1.2f $' % {"turb": tubulence},
       horizontalalignment='center',  verticalalignment='center',
       fontweight='bold', fontsize=TextFontSize)
p.grid(True)

p.subplot(313)
p.plot(data[:, 0], data[:, 2], 'k-', linewidth=1)
p.fill_between(data[:, 0], data[:, 2] + data[:, 3], data[:, 2] - data[:, 3],
               color=(0.75, 0.75, 0.75, 0.75))
p.axis([0, data[:, 0].max(), 0, (data[:, 2] + data[:, 3]).max()])
p.xlabel('time (steps)', fontsize=X_Label)
p.ylabel('mean number\n of genes', fontsize=Y_Label)
p.xticks(size=X_ticks)
p.yticks(size=Y_ticks)
p.grid(True)

p.show()
