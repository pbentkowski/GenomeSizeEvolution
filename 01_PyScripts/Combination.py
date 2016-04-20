#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the number of k-elenemt combinations in n-element set

Created on Wed Jan 11 18:05:48 2012

Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import numpy


def combinat(n, k):
    if n <= 0:
        print "Call from Combination.combinat(n, k) : n and k have to be",\
            "a positive intigers."
        return numpy.NaN
    elif k <= 0:
        print "Call from Combination.combinat(n, k) : n and k have to be",\
            "a positive intigers."
        return numpy.NaN
    elif n < k:
        print "Call from Combination.combinat(n, k) : n cannont be smaller",\
            "then k."
        return numpy.NaN
    elif (n % 1) != 0:
        print "Call from Combination.combinat(n, k) : n and k have to be",\
            "a positive intigers."
        return numpy.NaN
    elif (k % 1) != 0:
        print "Call from Combination.combinat(n, k) : n and k have to be",\
            "a positive intigers."
        return numpy.NaN
    else:
        return numpy.math.factorial(n) / (numpy.math.factorial(k)
                                          * numpy.math.factorial(n-k))
