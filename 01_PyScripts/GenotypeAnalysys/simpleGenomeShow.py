#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Tue Feb 26 17:41:09 2013
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""

import pylab as p
import re
import sys
import linecache as ln


def find_nearest(array, value):
    '''finds index of the nearest value'''
    return (p.absolute(array-value)).argmin()

FontSize = 22
TickSize = 19
Y_max = 1.0
Hpad = 0.2

GenerData = p.genfromtxt("GeneralData.dat")
try:
    theLine = find_nearest(GenerData[:, 0], int(sys.argv[1]))
    inpt = "Time step is " + str(int(sys.argv[1])) + " (line "\
        + str(theLine) + ")"
except:
    theLine = -1
    inpt = "Time step is " + str(int(GenerData[-1, 0])) + " (last line)"
genMeans = p.genfromtxt("GenotypeEnvelMean.dat")
genSTD = p.genfromtxt("GenotypeEnvelSTD.dat")
#GenerData = p.genfromtxt("GeneralData.dat")
l = re.split(" ", ln.getline("ModelParams.dat", 7))
xResolution = float(l[6])
l = re.split(" ", ln.getline("ModelParams.dat", 6))
T = float(l[6])
print inpt
xSize = 2.0 / xResolution
if(genMeans.shape[1] != xSize):
    print "ERROR: Env space size does not match genMeans space size! Check"
    exit()
x = p.arange(-1.0, 1.0, xResolution)

p.figure(1, figsize=(12, 8))
p.subplot2grid((3, 3), (0, 0), colspan=3)
try:
    p.vlines(GenerData[theLine, 0], -1.0, 1.0, 'k', linestyles='solid', lw=2)
except:
    p.vlines(GenerData[-1, 0], -1.0, 1.0, 'k', linestyles='solid', lw=2)
p.plot(GenerData[:, 0], GenerData[:, 1], 'k-')
p.axis([0, GenerData[:, 0].max(), -1, 1])
p.ylabel('environmental \n conditions ($ x $)', fontsize=FontSize)
p.xlabel('time (steps)', fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)

p.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=2)
p.plot(x, genMeans[theLine, :], linewidth=2, color='k')
p.fill_between(x, genMeans[theLine, :] + genSTD[theLine, :],
               genMeans[theLine, :] - genSTD[theLine, :],
               color=(0.75, 0.75, 0.75, 0.75))
p.ylim(0.0, Y_max)
p.xlabel('environmental conditions ($ x $)', fontsize=FontSize)
p.ylabel('avarege uptake efficiency ($ U(x) $)', fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)
p.tight_layout(h_pad=Hpad)
#p.tight_layout()
if theLine == -1:
    figFile = "/home/piotr/fig/" + str(T) + "_T.s_last.png"
else:
    figFile = "/home/piotr/fig/" + str(T) + "_T.s" + str(theLine) + ".png"
p.savefig(figFile, dpi=150)
p.show()
