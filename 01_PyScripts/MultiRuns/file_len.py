#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Just an utility function. Counts the number of lines in a text file. 

Created on Mon May 30 17:38:33 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1