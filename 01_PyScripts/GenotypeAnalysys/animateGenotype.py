#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Tue Feb 26 22:13:39 2013
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""
import pylab as p
import re
import linecache as ln
import matplotlib.animation as animation

everyOtherRow = 10

FontSize = 20
TickSize = 17
Y_max = 1.0
Hpad = 0.05
p.rc('xtick', labelsize=TickSize)
p.rc('ytick', labelsize=TickSize)


genMeans = p.genfromtxt("GenotypeEnvelMean.dat")
genSTD = p.genfromtxt("GenotypeEnvelSTD.dat")
GenerData = p.genfromtxt("GeneralData.dat")
l = re.split(" ", ln.getline("ModelParams.dat", 7))
xResolution = float(l[6])
l = re.split(" ", ln.getline("ModelParams.dat", 6))
T = float(l[6])
print "Turbulence level = ", T
xSize = 2.0 / xResolution
if(genMeans.shape[1] != xSize):
    print "ERROR: Env space size does not match genMeans space",
    print "size! Check them."
    exit()
x = p.arange(-1.0, 1.0, xResolution)

# -- trimmig rows
fullPlot = GenerData[:, 0:2]
lastRow = genMeans[-1, :]
genMeans = genMeans[::everyOtherRow]
genMeans = p.vstack([genMeans, lastRow])
lastRow = genSTD[-1, :]
genSTD = genSTD[::everyOtherRow]
genSTD = p.vstack([genSTD, lastRow])
lastRow = GenerData[-1, :]
GenerData = GenerData[::everyOtherRow]
GenerData = p.vstack([GenerData, lastRow])
# -- trimmed

fig = p.figure(1, figsize=(15, 9))
ax1 = p.subplot2grid((3, 3), (0, 0), colspan=3)
ax1.axis([0, fullPlot[:, 0].max(), -1, 1])
ax1.set_ylabel('environmental \n conditions ($ x $)', fontsize=FontSize)
ax1.set_xlabel('time (steps)', fontsize=FontSize)
ax1.plot(fullPlot[:, 0], fullPlot[:, 1], 'k-')
ax1.grid(True)

ax2 = p.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=2)
#ax2.axis([-1.0, 1.0, 0, 1.0])
#ax2.set_xlabel('environmental conditions ($ x $)', fontsize=FontSize)
#ax2.set_ylabel('average uptake efficiency ($ U(x) $)', fontsize=FontSize)
#ax2.grid(True)
p.tight_layout(h_pad=Hpad)

line1, = ax1.plot([], [], 'ro', ms=10)
line3, = ax2.plot([], [], linewidth=2, color='k')


def init():
    line1.set_data([])
    return line1


def animate(i):
    line1.set_data(GenerData[i, 0], GenerData[i, 1])
    ax2.cla()
    ax2.axis([-1.0, 1.0, 0, 1.0])
    ax2.set_xlabel('environmental conditions ($ x $)', fontsize=FontSize)
    ax2.set_ylabel('average uptake efficiency ($ U(x) $)', fontsize=FontSize)
    p.fill_between(x, genMeans[i, :] + genSTD[i, :],
                   genMeans[i, :] - genSTD[i, :],
                   color=(0.75, 0.75, 0.75, 0.75))
    p.plot(x, genMeans[i, :], linewidth=2, color='k')
    p.grid(True)
    p.tight_layout(h_pad=Hpad)
    return line1

ani = animation.FuncAnimation(fig, animate, interval=0.01, blit=False)
#ani.save('genome_evol.mp4', fps=25)
p.show()
