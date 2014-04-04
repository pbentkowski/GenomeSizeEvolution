#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 18:37:02 2010

@author: Piotr Bentkowski
"""

import pylab as p

def gaussar(x, params):
        """ gaussar(x,[amplitude, sigma, c]) Returns gaussian curve over x with 
	maximum value of amplitude at c and with given sigma"""
        y = params[0] * p.exp( -(x-params[2])**2/(2*params[1]**2 ) )
        return y
        
