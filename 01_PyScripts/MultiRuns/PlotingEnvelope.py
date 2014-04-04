#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 12:18:54 2010

@author: Piotr Bentkowski
"""
import pylab as p
import matplotlib.cm as cm

data = p.genfromtxt("GeneralData.dat")
envelopes_sizes = p.genfromtxt("FrameSizeData.dat")
envelopes_maximums = p.genfromtxt("FrameMaxData.dat")
rel_envelopes_sizes = p.zeros((envelopes_sizes.shape[0],
                              envelopes_sizes.shape[1]))
rel_envelopes_maximums = p.zeros((envelopes_maximums.shape[0],
                                  envelopes_maximums.shape[1]))

# calculating stuff for envelope size distribution plot
for i in range(envelopes_sizes.shape[0]):
        S = sum(envelopes_sizes[i, :])
        for j in range(envelopes_sizes.shape[1]):
            rel_envelopes_sizes[i, j] = envelopes_sizes[i, j] / S

# calculating stuff for maximums of distribution plot
for k in range(envelopes_maximums.shape[0]):
        MaxSum = sum(envelopes_maximums[k, :])
        for l in range(envelopes_maximums.shape[1]):
            rel_envelopes_maximums[k, l] = envelopes_maximums[k, l] / MaxSum

#---- sixth plot-------
step = 1.0/envelopes_sizes.shape[1]
y = p.arange(0.0, 1.0, step)
p.figure(6, figsize=(12, 6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.contourf(data[:, 0], y, rel_envelopes_sizes.transpose(), cmap=cm.gist_yarg)
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)')
p.ylabel('envelopes\' surfeace / evironmetal space')
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence')
#---- end of sixth plot-------

#---- seventh plot-------
step_a = 1.0/envelopes_maximums.shape[1]
yy = p.arange(0.0, 1.0, step_a)
p.figure(7, figsize=(12, 6))
#p.get_current_fig_manager().window.wm_geometry("1000x600+300+160")
p.contourf(data[:, 0], yy, rel_envelopes_maximums.transpose(),
           cmap=cm.gist_yarg)
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  y.max()])
p.xlabel('time (steps)')
p.ylabel('envelopes\' maximums')
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence')
#---- end of seventh plot-------

p.show()
print "done!"
