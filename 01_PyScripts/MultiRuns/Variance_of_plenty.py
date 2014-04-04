#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter the doc string...

Created on Tue Feb 28 01:09:09 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import os
import pylab as p

AxisLabelFontSize = 22
AxisTickFontSize = 22
AnnotateFontSize = 19

dec_places = '%1.4f'

def LoadEnvData(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname,'VarianceOfEnvEtAl.dat'):
            data = p.genfromtxt(filepath)
            # goes like this: mean envelope surface, turbulance level T,
            # STD of envelope surface
            arg.append((data[:,2], data[:,3], data[:,0]))

EnvData = []
os.path.walk(os.getcwd(), LoadEnvData, EnvData)
print EnvData[0]

p.figure(1, figsize=(1000, 600))
i = 0
for item in EnvData:
    ax = p.plot(item[0], item[1], 'k-')
    i = i + 1
p.ylabel('variance of env contitions ($ s^{2}(x) $)', fontsize = AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize = AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid()
p.show()