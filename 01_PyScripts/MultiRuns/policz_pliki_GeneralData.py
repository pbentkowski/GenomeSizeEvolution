#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 15:52:07 2010

@author: piotr
"""

import os

def CountFiles(arg, dirname, files):
    counter = 0
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname,'GeneralData.dat'):
            counter += 1
    arg.append((counter))

Data = []
os.path.walk(os.getcwd(), CountFiles, Data)

print "There are", sum(Data), "model run results in this directory."
