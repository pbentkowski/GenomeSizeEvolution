#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Enter the doc string...

Created on Mon Mar  5 01:38:54 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import pylab as p

LabelFontSize=20
TickSize=15

TheMeans = p.genfromtxt("VariancesFrame.dat")

p.figure(9, figsize=(1000,600))
p.plot(TheMeans[:, 0], TheMeans[:, 2], 'ko-')
p.fill_between(TheMeans[:, 0], TheMeans[:, 2], 0.0 ,
               color=(0.75, 0.75, 0.75, 0.75))
p.xlabel('turbulence level ($ T $)', fontsize=LabelFontSize)
p.ylabel("variance of the env conditions",
         fontsize=LabelFontSize)
p.ylim(ymin=0)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.axis([0, 0.5, 0, 0.60001])
p.grid(True)
p.show()
