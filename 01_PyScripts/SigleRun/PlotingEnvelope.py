#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 12:18:54 2010

@author: Piotr Bentkowski
"""
import pylab as p
import matplotlib.cm as cm
import re
import linecache as ln

LabelsFontSize = 22
AxisTickFontSize = 20

par_0 = re.split(" ", ln.getline('ModelParams.dat', 6))
tubulence = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 32))
numb_of_bins_envel = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 33))
bin_width_envel = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 34))
numb_of_bins_max = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 35))
bin_width_max = float(par_0[6])

data = p.genfromtxt("GeneralData.dat")
envelopes_sizes = p.genfromtxt("FrameSizeData.dat")
envelopes_maximums = p.genfromtxt("FrameMaxData.dat")
rel_envelopes_sizes= p.zeros((envelopes_sizes.shape[0],envelopes_sizes.shape[1]))
rel_envelopes_maximums = p.zeros((envelopes_maximums.shape[0],
                                  envelopes_maximums.shape[1]))

# calculating stuff for envelope size distribution plot
for i in range(envelopes_sizes.shape[0]):
        S = sum(envelopes_sizes[i, :])
        for j in range(envelopes_sizes.shape[1]):
            rel_envelopes_sizes[i,j] = envelopes_sizes[i, j] / S

# calculating stuff for maximums of distribution plot
for k in range(envelopes_maximums.shape[0]):
        MaxSum = sum(envelopes_maximums[k, :])
        for l in range(envelopes_maximums.shape[1]):
            rel_envelopes_maximums[k,l] = envelopes_maximums[k,l] / MaxSum

#---- sixth plot-------
y = p.arange(0.0, numb_of_bins_envel * bin_width_envel, bin_width_envel)
p.figure(6, figsize=(1000,600))
p.contourf(data[:, 0], y, rel_envelopes_sizes.transpose(), cmap=cm.gist_yarg)
p.axis([data[:, 0].min(), data[:, 0].max(), y.min(), y.max()+bin_width_envel])
p.xlabel('time (steps)', fontsize=LabelsFontSize)
p.ylabel('envelope\'s surfeace to \n evironmetal space ratio',
         fontsize=LabelsFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence', fontsize=LabelsFontSize)
#---- end of sixth plot-------

#---- seventh plot-------
yy = p.arange(0.0, numb_of_bins_max * bin_width_max, bin_width_max)
p.figure(7, figsize=(1000,600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+300+160")
p.contourf(data[:, 0], yy, rel_envelopes_maximums.transpose(), cmap=cm.gist_yarg)
p.axis([data[:, 0].min(),  data[:, 0].max(),  y.min(),  1.0])
p.xlabel('time (steps)', fontsize=LabelsFontSize)
p.ylabel('envelopes\' maximums', fontsize=LabelsFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
cb = p.colorbar(format=r"%.2f")
cb.set_label('frequency of occurrence', fontsize=LabelsFontSize)
p.title('Turbulence level = %(turb)1.3f'% {"turb":tubulence},
        fontsize=LabelsFontSize)

#---- end of seventh plot-------
p.show()
print "done!"