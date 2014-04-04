#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Return amount of digits of given number. For complex numbers and strings
returnsNumPy.NaN

Created on Tue Jan 31 18:15:41 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import numpy

      
def digits(x):
    """ Returns amount of digits of x. """
    try:
        float(x)
    except TypeError:
        return numpy.nan
    except ValueError:
        return numpy.nan
    x = numpy.absolute(x)
    return int(numpy.math.ceil(numpy.math.log10(x)))
