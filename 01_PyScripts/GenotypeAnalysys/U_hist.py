#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter a doc string here.....

------------------------------------------------------
Created on Fri Mar  1 15:05:11 2013
@author: Piotr Bentkowski - bentkowski.piotr@gmail.com
------------------------------------------------------
"""

import pylab as p
import re
import linecache as ln

theLine = -1

hist = p.genfromtxt("AvaregeGenotype.dat")
l = re.split(" ", ln.getline("ModelParams.dat", 7))
U_res  = float(l[6])
numOfBins = 1.0 / U_res
x = p.arange(0.0, 1.0, U_res)

hist[:,0] = 0.0
p.figure(1, figsize=(10, 6))
p.bar(x, hist[-1,:], width=U_res, color='k', linewidth=0)
p.grid(True)
p.show()