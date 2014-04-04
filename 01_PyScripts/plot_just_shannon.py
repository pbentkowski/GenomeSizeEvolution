#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Put a doc string here...

Created on Thu Mar  8 17:58:02 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import pylab as p

AxisLabelFontSize = 22
AxisTickFontSize = 22

data = p.genfromtxt("GeneralData.dat")
p.figure(num=None, figsize=(1000,600))
p.subplot(211)
p.plot(data[:, 0], data[:, 9], 'k-')
p.axis([data[:, 0].min(),  data[:, 0].max(), 0.9*data[:, 9].min(),
        1.1*data[:, 9].max()])
p.xlabel('time (steps)', fontsize=AxisLabelFontSize)
p.ylabel('Shannon index', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
p.show()