#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Put a doc string here...

Created on Fri Mar  2 18:02:23 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
#import pylab as p
import os
import re
import linecache as ln


def ReplaceOldTurb(arg, dirname, files):
        for file in files:
            filepath = os.path.join(dirname, file)
            if filepath == os.path.join(dirname, 'ModelParams.dat'):
                l = re.split(" ", ln.getline(filepath, 6))
                turb_param = float(l[6])
                numm = turb_param * 0.5
                new_string = '    turbulence = %(num)1.4f\n' % {"num": numm}
                lines = open(filepath, 'r').readlines()
                lines[5] = new_string
                out = open(filepath, 'w')
                out.writelines(lines)
                out.close()
                arg.append((turb_param, numm))

CheckData = []
os.path.walk(os.getcwd(), ReplaceOldTurb, CheckData)

print "Done!"
