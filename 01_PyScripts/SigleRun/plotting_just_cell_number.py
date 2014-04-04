#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Put a doc string here...

Created on Sat Sep 29 20:28:46 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import pylab as p
import re
import linecache as ln
#import matplotlib.cm as cm

the_cool_line = 53
data = p.genfromtxt("GeneralData.dat")
l_1 = re.split(" ", ln.getline('ModelParams.dat', the_cool_line))
The_big_value = float(l_1[6])
l_1 = re.split(" ", ln.getline('ModelParams.dat', 6))
print "T_value =", float(l_1[6])

#--- Fonts sizes ---
X_ticks = 15
Y_ticks = 15
X_Label = 18
Y_Label = 18
AnnotateFontSize = 18
XX = 180000
YY = 500
#---- first plot-------
p.figure(1, figsize=(10, 6))

#p.subplot(311)
#p.plot(data[:, 0], data[:, 1],'k-')
#p.axis([0, data[:,0].max(), -1, 1])
#p.ylabel('environmental conditions', fontsize=X_Label)
#p.xticks(size=X_ticks)
#p.yticks(size=Y_ticks)
#p.grid(True)

p.subplot(313)
p.plot(data[:, 0],data[:, 4],'k-')
p.annotate('$h_{c} = %1.4f$'%(The_big_value,), xy = (XX, YY),  xycoords='data',
           fontsize=AnnotateFontSize)
p.axis([0, data[:,0].max(), 0, data[:,4].max()])
p.xlabel('time (steps)', fontsize=X_Label)
p.ylabel('number of cells', fontsize=X_Label)
p.xticks(size=X_ticks)
p.yticks(size=Y_ticks)
p.grid(True)

#p.subplot(313)
#p.plot(data[:, 0], data[:, 2],'k-', linewidth=1)
#p.fill_between(data[:, 0],data[:, 2] + data[:, 3],data[:, 2] - data[:, 3],
#               color=(0.75,0.75,0.75,0.75))
#p.axis([0, data[:,0].max(), 0, (data[:, 2] + data[:, 3]).max()])
#p.xlabel('time (steps)', fontsize=X_Label)
#p.ylabel('mean number of genes', fontsize=Y_Label)
#p.xticks(size=X_ticks)
#p.yticks(size=Y_ticks)
#p.grid(True)
#---- end of first plot-------
p.show()