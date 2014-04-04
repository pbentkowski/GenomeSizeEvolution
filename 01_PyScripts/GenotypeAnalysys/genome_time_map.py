#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Wed Feb  5 01:27:15 2014
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""
import re
import pylab as p
import linecache as ln

FontSize = 32
TickSize = 28
Y_max = 1.0
Hpad = 0.2

GenerData = p.genfromtxt("GeneralData.dat")
data = p.genfromtxt("GenotypeEnvelMean.dat")
l = re.split(" ", ln.getline("ModelParams.dat", 7))
xResolution = float(l[6])
l = re.split(" ", ln.getline("ModelParams.dat", 6))
T = float(l[6])
print "Turbulence level =", T
dM = p.ma.masked_greater(data, 1.0)
x = GenerData[:, 0]
y = p.linspace(-1.0, 1.0, data.shape[1])

p.figure(1, figsize=(30, 6))
p.subplot(121)
p.plot(GenerData[:, 0], GenerData[:, 1], 'k-', lw=2)
p.axis([0, GenerData[:, 0].max(), -1, 1])
#p.ylabel('environmental \n conditions ($ x $)', fontsize=FontSize)
#p.xlabel('time (steps)', fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)

p.subplot(122)
p.pcolormesh(x, y, dM.transpose(), cmap=p.cm.gray_r, vmin=0, vmax=1)
#p.xlabel('time (steps)', fontsize=FontSize)
#p.ylabel('   environmental \n conditions ($ x $)', fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)
cb = p.colorbar(format=r"%.1f", orientation='horizontal')
cb.ax.tick_params(labelsize=TickSize)
cb.set_label(r'uptake efficiency $U(x)$', fontsize=FontSize)
cb.ax.minorticks_on()
p.tight_layout(h_pad=Hpad)
figFile = "/home/piotr/fig/" + str(T) + "_T_time_map.png"
p.savefig(figFile, dpi=150)
p.show()
