#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 23:40:03 2010

@author: piotr
"""

import os, re
import pylab as p
import linecache as ln

def LoadGeneNumbers(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname,'GenomeSizes.dat'):
            genes = p.genfromtxt(filepath)
            arg.append((genes[:,0], genes[:,2]))
            
def LoadGeneCost(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname,'Evelopes.dat'):
            l = re.split(" ", ln.getline(filepath, 1))
            gene_cost = float(l[1])              
            arg.append((gene_cost))         
            
               

AllData = []
os.path.walk(os.getcwd(), LoadGeneNumbers, AllData)
Tags = []
os.path.walk(os.getcwd(), LoadGeneCost, Tags)
i = 0 
p.figure(1, figsize=(1000,600))
for item in AllData:
    XX = float(item[1][-i-2])
    YY = float(item[0][33])
    ax = p.plot(item[1], item[0], 'k-')
    ax = p.annotate('%4.2f'%(Tags[i],), xy = (XX, YY),  xycoords='data', fontsize=14)
    i = i + 1
p.xlabel('turbulence level ($ T $)', fontsize=17)
p.ylabel('mean number of genes', fontsize=17)
p.xticks(size=16)
p.yticks(size=16)
p.grid()
p.show()