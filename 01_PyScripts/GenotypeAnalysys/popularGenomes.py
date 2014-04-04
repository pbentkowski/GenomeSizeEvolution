#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Thu Feb 28 21:33:57 2013
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""
import pylab as p
import re
import linecache as ln
import gaussar as ga  # this one is mine

FontSize = 24
TickSize = 21

genMeans = p.genfromtxt("GenotypeEnvelMean.dat")
l = re.split(" ", ln.getline("ModelParams.dat", 7))
xResolution = float(l[6])
l = re.split(" ", ln.getline("ModelParams.dat", 6))
T = float(l[6])
xSize = 2.0 / xResolution
if(genMeans.shape[1] != xSize):
    print "ERROR: Env space size does not match genMeans space size! Check"
    exit()

with open("GenomesOfPopulation.dat") as f:
    theData = f.read().splitlines()

counter = 0
line_numb = 0
NumOfCells = 0.0
for line in theData:
    try:
        l = re.split(" ", line)
        number = float(l[0])
        NumOfCells = NumOfCells + number
        if(counter <= number):
            counter = number
            crutial_line = line_numb
    except:
        pass
    line_numb = line_numb + 1
#print crutial_line
print "There are", int(NumOfCells), "cells in the population"

L2 = re.split(" ", theData[crutial_line])
rows = (len(L2) - 4) / 3
genes = p.zeros((rows, 3))
if rows is 1:
    descr = "T = " + str(T) + "; " + str(rows) + " genes; "
    perct = counter/NumOfCells * 100
    PERC = '%(perc)1.2f %% of population' % {"perc": perct}
    descr = descr + PERC
    print descr
else:
    descr = "T = " + str(T) + "; " + str(rows) + " genes; "
    perct = counter/NumOfCells * 100
    PERC = "%(perc)1.2f %% of population" % {"perc": perct}
    descr = descr + PERC
    print descr
print "Hight (A) , sigma , positions of max (c)"
for i in xrange(len(L2) - 4):
    ii = p.floor(i / 3)
    jj = i % 3
    if jj is 0:
        genes[ii, 0] = float(L2[i+3])
    elif jj is 1:
        genes[ii, 1] = float(L2[i+3])
    elif jj is 2:
        genes[ii, 2] = float(L2[i+3])
    else:
        print "ERROR in loading genes to numpy array."
        exit()
print genes

x = p.arange(-1.0, 1.0, 0.01)
xx = p.arange(-1.0, 1.0, xResolution)
p.figure(1, figsize=(13, 7))
p.axis([-1, 1, 0, 1])
p.xlabel('environmental conditions ($ x $)', fontsize=FontSize)
p.ylabel('uptake efficiency ($ U(x) $)', fontsize=FontSize)
p.title(descr, fontsize=FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)
p.plot(xx, genMeans[-1, :], 'k-', lw=4)  # , color=(0.75, 0.75, 0.75, 1.0))
for k in xrange(genes.shape[0]):
    gns = ga.gaussar(x, genes[k, :])
#    p.plot(x, gns, 'k')
    p.fill_between(x, gns, facecolor=(0.75, 0.75, 0.75, 0.75))
out_name = "/home/piotr/fig/cool_genome_" + str(T) + ".png"
p.savefig(out_name, dpi=150)
p.show()
