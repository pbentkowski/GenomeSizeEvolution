#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 16:33:22 2010

@author: Piotr Bentkowski
"""
import pylab as p

#Mean_Var = p.genfromtxt("ENV_variance.dat")
#tags = p.genfromtxt("ENV_TAGS_variance.dat")
#All_the_frames = p.genfromtxt("ENV_All_the_frames.dat")


tags = p.arange(0.0,1.0,0.01)
half_tags = tags * 0.5
All_the_frames = p.arange(2,200,2)
ENV = 2.0 * p.ones((25000, tags.shape[0]))

for i in xrange(ENV.shape[1]):
    for j in xrange(1,ENV.shape[0]):
        ENV[0,i] = 0.0
        while ENV[j,i] > 1.0 or ENV[j,i] < -1.0:
              ENV[j,i] = ENV[j-1,i] + tags[i]*p.random() - half_tags[i]

# --- caunting the variances ---
Mean_Var = p.zeros((All_the_frames.shape[0] , ENV.shape[1]))
for i in xrange(ENV.shape[1]):
        print 'Processing data series No.',i+1, 'out of', ENV.shape[1]
        l = 0
        for item in All_the_frames:
            diff_of_means_size = ENV.shape[0]- item
            diff_of_means = p.zeros((diff_of_means_size, 1))
            half = item/2
            k = half
            for j in range(diff_of_means.shape[0]):
                diff_of_means[j] = ENV[k-half:k+half, i].var()
                k = k + 1
            Mean_Var[l,i] = diff_of_means.sum() / float(diff_of_means_size)
            l = l + 1
        i = i + 1
print 'I\'m done with processing data series. Writing to files.'

p.savetxt("ENV_variance.dat", Mean_Var)
p.savetxt("ENV_TAGS_variance.dat", tags)
p.savetxt("ENV_All_the_frames.dat", All_the_frames)
