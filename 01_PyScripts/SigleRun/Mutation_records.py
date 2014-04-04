#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Opens the file MutationRecords.dat and calculates the mean and the STD of numbers
of mutations of all three kinds separately. Takes under account that there are
different clones and the number of mutations is averaged over a number of clones.

Created on Tue Jul 19 12:27:19 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import linecache as ln
import re
# the one below is mine
import file_len as fl

def MutationPerClone(x):
    '''Calculates average number and STD of mutations from MutationRecords.dat
    file per clone. Mut_record[0] - number of cells, [1] - number of clones, [2]
    - mean number of point mutations, [3] - STD of point mutations, [4] - mean
    number of duplications, [5] - STD  of duplications, [6] - mean number of
    deletion, [7] - STD of deletion.'''
    Mut_record = p.zeros(8)
    Mut_record[0] = x.shape[0]
    ZERRO = p.zeros((x.shape[0], 1))
    x = p.concatenate((x, ZERRO), axis=1)
    for i in xrange(0, x.shape[0]):
        if (x[i, 4] == 0):
            for j in xrange(i+1, x.shape[0]):
                if (x[j, 0] == x[i, 0]):
                    x[j, 4] = p.nan
    x = x[~p.isnan(x).any(1)]
    Mut_record[1] = x.shape[0]
    Mut_record[2] = x[:,1].mean()
    Mut_record[3] = x[:,1].std()
    Mut_record[4] = x[:,2].mean()
    Mut_record[5] = x[:,2].std()
    Mut_record[6] = x[:,3].mean()
    Mut_record[7] = x[:,3].std()
    return Mut_record

l = re.split(" ", ln.getline('MutationRecords.dat',
                             fl.file_len('MutationRecords.dat')-2))
is_alive = float(l[0])
if (is_alive != -1):
    # Printing stuff just in case...
    mutations = p.genfromtxt('MutationRecords.dat')
    MUTT = MutationPerClone(mutations)
    print "Number of cells :", MUTT[0]
    print "Number of clonal strains :", MUTT[1]
    print "---------------------------------------------"
    print "Average number of gene modification:", MUTT[2], "+/-", MUTT[3], \
    "(mean +/- STD)"
    print "Average number of duplications     :", MUTT[4], "+/-", MUTT[5], \
    "(mean +/- STD)"
    print "Average number of deletions        :", MUTT[6], "+/-", MUTT[7], \
    "(mean +/- STD)"
else:
     print "This run didn\'t make it to the end."
     exit()