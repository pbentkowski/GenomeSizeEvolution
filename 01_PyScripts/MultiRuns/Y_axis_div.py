#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Divides Y axis into a given number of ticks.

Created on Wed Feb  1 01:03:24 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""

import numpy as p
import Digits as d # this one is mine

def divide_Y_axis(MIN, MAX, NUM_OF_TICKS):
    """Divides Y axis into a given number of ticks. MIN is the absolute
    beginning. MAX is the maximal value of data set you have. The true maximum
    of the axes gests adjusted by this function."""
    if MAX <= MIN:
        print "Error call from Y_axis_div.divide_Y_axis: MAX has to be"\
        " larger then MIN."
        return [p.nan, p.nan, p.nan]
    elif NUM_OF_TICKS <= 1:
        print "Error call from Y_axis_div.divide_Y_axis: NUM_OF_TICKS"\
        " has to be larger then 1."
    else:
        NUM_OF_TICKS = p.round(NUM_OF_TICKS)
        MAX = MAX - MIN
        scaling_fac = 10.0**(d.digits(MAX) - 1)
        if (MAX / scaling_fac) < (NUM_OF_TICKS - 1.0):
           scaling_fac = 10.0**(d.digits(MAX) - 2)    
        relative_max = p.ceil(MAX / scaling_fac) * scaling_fac
        interv = p.ceil(relative_max / (NUM_OF_TICKS * 10.0 \
            ** (d.digits(scaling_fac)))) * 10.0 ** (d.digits(scaling_fac))
        dig_one = d.digits(interv)
        dig_two = d.digits(2.0 * interv)
        if dig_one < dig_two:
            if (interv % (5.0 * (dig_two - 1))):
                interv = 10.0 ** (d.digits(2.0 * interv) - 1)
        relative_max = NUM_OF_TICKS * interv
        return [MIN, relative_max, interv]