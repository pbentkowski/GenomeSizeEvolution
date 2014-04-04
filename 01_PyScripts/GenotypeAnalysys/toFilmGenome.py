#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Sat Mar  2 17:47:33 2013
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""
import pylab as p
import re
import linecache as ln


everyOtherRow = 2

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
    print "ERROR: Env space size does not match genMeans space size!",
    print "Check them."
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
ax2 = p.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=2)

numbOfFrames = genMeans.shape[0]
for i in xrange(numbOfFrames):
    if i < 10:
        nameFile = "0000" + str(i)
    elif i >= 10 and i < 100:
        nameFile = "000" + str(i)
    elif i >= 100 and i < 1000:
        nameFile = "00" + str(i)
    elif i >= 1000 and i < 10000:
        nameFile = "0" + str(i)
    else:
        nameFile = str(i)
    ax1.cla()
    ax2.cla()
    ax1 = p.subplot2grid((3, 3), (0, 0), colspan=3)
    ax1.axis([0, fullPlot[:, 0].max(), -1, 1])
    ax1.set_ylabel('environmental \n conditions ($ x $)', fontsize=FontSize)
    ax1.set_xlabel('time (steps)', fontsize=FontSize)
    ax1.grid(True)
#    ax1.vlines(GenerData[i, 0], -1.0, 1.0, 'r', lw=2)
    ax1.plot(fullPlot[:, 0], fullPlot[:, 1], 'k-')
    ax1.plot(GenerData[i, 0], GenerData[i, 1], 'ro', ms=10)
    ax1.vlines(GenerData[i, 0], -1.0, 1.0, 'r', lw=2)
    ax2 = p.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=2)
    ax2.axis([-1.0, 1.0, 0, 1.0])
    ax2.set_xlabel('environmental conditions ($ x $)', fontsize=FontSize)
    ax2.set_ylabel('average uptake efficiency ($ U(x) $)', fontsize=FontSize)
    alpha = 1.0
    if i > 9:
        for j in xrange(10):
            ax2.vlines(GenerData[i-j, 1], 0.0, 1.0,
                       lw=2, color=(1.0, 0, 0, alpha))
            alpha = alpha * 0.7
    else:
        ax2.vlines(GenerData[i, 1], 0.0, 1.0,
                   lw=2, color=(1.0, 0, 0, alpha))
    ax2.fill_between(x, genMeans[i, :] + genSTD[i, :],
                     genMeans[i, :] - genSTD[i, :],
                     color=(0.75, 0.75, 0.75, 0.75))
    ax2.plot(x, genMeans[i, :], linewidth=2, color='k')
    ax2.grid(True)
    fig.tight_layout(h_pad=Hpad)
    nameFile = nameFile + ".png"
    fig.savefig(nameFile)
    print "Done with frame", i, "out of", numbOfFrames - 1
print "All frames done! Processing movie."

# Generating video out of it
#ffmpeg -i %05d.png -f avi -vcodec mjpeg -qscale 1 output.avi
bash = "ffmpeg -i %05d.png -f avi -vcodec mjpeg -vf scale=1080:-1 output.avi"
import subprocess
process = subprocess.Popen(bash.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]
print "All done!"
