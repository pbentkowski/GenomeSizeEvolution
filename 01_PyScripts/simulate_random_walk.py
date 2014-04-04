#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simulates the random walk used to generate variable environment in real 1:1
resolution. Remember, that in the model's implementation only every 10th step
is saved to file.

Created on Wed Feb 29 00:38:50 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""

import pylab as p

N = 200000
env = 0.0
turbulence = p.array([0.000, 0.005, 0.007, 0.010, 0.012, 0.015, 0.017, 0.020,\
 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.200, 0.300, 0.400,\
 0.500, 0.600, 0.700, 0.800, 0.900, 1.000])

the_array = p.zeros((turbulence.shape[0],N))
varrrr = p.zeros(turbulence.shape[0])

for j in xrange(turbulence.shape[0]):
    half_of_turbulence = turbulence[j] * 0.5
    for i in range(N):
        env = env + (turbulence[j] * p.random() - half_of_turbulence)
        while (env * env) > 1.0:
            env = env + (turbulence[j] * p.random() - half_of_turbulence)
        the_array[j, i] = env
    print 'done for turbulence:', turbulence[j]
    varrrr[j] = the_array[j, :].var()


print the_array.shape[0]

p.figure(num=None, figsize=(8, 6), dpi=80)
p.plot(varrrr)
p.grid(True)
p.show()