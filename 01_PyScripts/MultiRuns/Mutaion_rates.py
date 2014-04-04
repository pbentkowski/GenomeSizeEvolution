#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Iterates over a bunch o directories and opens the files MutationRecords.dat.
Calculates the mean and the STD of numbers of mutations of all three kinds
separately. Takes under account that there are different clones and the number
of mutations is averaged over a number of clones. Writes the results to a file
MutationStats.dat.

Created on Tue Jul 12 14:27:53 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import os, re
import linecache as ln
# the one below is mine
import file_len as fl

AxisLabelFontSize = 17
AxisTickFontSize = 17
AnnotateFontSize = 15

MaxMutationNumber = 70.0
Letter_location = MaxMutationNumber - 5.0

time_frame = 1e5

def MutationPerClone(x):
    '''Calculates average number and STD of mutations from MutationRecords.dat
    file per clone. Mut_record[0] - number of cells, [1] - number of clones, [2]
    - mean number of point mutations, [3] - STD of point mutations, [4] - mean
    number of duplications, [5] - STD  of duplications, [6] - mean number of
    deletion, [7] - STD of deletion.'''
    Mut_record = p.zeros(8)
    Mut_record[0] = x.shape[0]
    ZERRO = p.zeros((x.shape[0],1))
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


def LoadMutationData(arg, dirname, files):
    '''Walks thorough directory tree and looks into MutationRecords.dat files
    calculating the average number and STD of mutations per clone, but only if
    population makes it to the end of simulation. Append all thet that to a
    list'''
    for file in files:
        filepath_general = os.path.join(dirname, file)
        if filepath_general == os.path.join(dirname,'GeneralData.dat'):
            l = re.split(" ", ln.getline(filepath_general,
                                         fl.file_len(filepath_general)-1))
            is_alive = float(l[4])
            if  is_alive > 0.0:
                filepath_turb = os.path.join(dirname,'ModelParams.dat')
                l = re.split(" ", ln.getline(filepath_turb, 6))
                turb_param = float(l[6])
                # change the line number in ln.getline below to extract the
                # parameter you're interested in
                l = re.split(" ", ln.getline(filepath_turb, 12))
                interesting_param = float(l[6])
                filepath_mut = os.path.join(dirname,'MutationRecords.dat')
                Mut_raw_data = p.genfromtxt(filepath_mut)
                l_2 = re.split(" ", ln.getline(filepath_mut,
                                               fl.file_len(filepath_mut)-1))
                time_stamp = float(l_2[1])
                arg.append((turb_param, MutationPerClone(Mut_raw_data),
                            interesting_param, time_stamp))
            else:
                break

MutationData = []
os.path.walk(os.getcwd(), LoadMutationData, MutationData)
MutationData = sorted(MutationData, key=lambda data_sort: data_sort[0])

TheFinalArray = p.zeros((len(MutationData),12))
i = 0
for item in MutationData:
    scale_fac = item[3] / time_frame
    TheFinalArray[i, 0] = item[0]    # Turbulence parameter
    TheFinalArray[i, 1] = item[1][0] # number of cells
    TheFinalArray[i, 2] = item[1][1] # number of clonal strains
    TheFinalArray[i, 3] = item[1][2] / scale_fac # mean number of point mutations
    TheFinalArray[i, 4] = item[1][3] / scale_fac # STD of number of point mut.
    TheFinalArray[i, 5] = item[1][4] / scale_fac # mean number of duplications
    TheFinalArray[i, 6] = item[1][5] / scale_fac # STD of number of duplications
    TheFinalArray[i, 7] = item[1][6] / scale_fac # mean number of deletion
    TheFinalArray[i, 8] = item[1][7] / scale_fac # STD of number of deletion
    TheFinalArray[i, 9] = item[2]    # interesting param (usually mutarion rate)
    TheFinalArray[i, 10] = item[3]   # after how many time steps was calculated
    TheFinalArray[i, 11] = time_frame # reference time frame
    i += 1

p.savetxt("MutationStats.dat", TheFinalArray, fmt='%1.5f')

p.figure(1, figsize=(1000, 600))
tick_marks_01 = p.arange(0., MaxMutationNumber + 5., 10.)
p.subplot(221)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 3], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 3] - TheFinalArray[:, 4],
               TheFinalArray[:, 3] + TheFinalArray[:, 4],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.95, Letter_location, 'A', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('Number of gene modifications', fontsize = AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_01, size=AxisTickFontSize)
p.ylim(ymax=MaxMutationNumber)
p.xlim(xmin=0.0)
p.grid(True)

p.subplot(222)
tick_marks_02 = p.arange(0., MaxMutationNumber + 5., 10.)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 5], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 5] - TheFinalArray[:, 6],
               TheFinalArray[:, 5] + TheFinalArray[:, 6],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.95, Letter_location, 'B', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('Number of duplications', fontsize = AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_02, size=AxisTickFontSize)
p.ylim(ymax=MaxMutationNumber)
p.xlim(xmin=0.0)
p.grid(True)

p.subplot(223)
tick_marks_03 = p.arange(0., MaxMutationNumber + 5., 10.)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 7], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 7] - TheFinalArray[:, 8],
               TheFinalArray[:, 7] + TheFinalArray[:,8],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.95, Letter_location, 'C', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('Number of deletions', fontsize = AxisLabelFontSize)
p.xlabel('Turbulence level ($ T $)', fontsize = AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_03, size=AxisTickFontSize)
p.ylim(ymax=MaxMutationNumber)
p.xlim(xmin=0.0)
p.grid(True)

p.subplot(224)
tick_marks_04 = p.arange(0., p.amax(TheFinalArray[:, 2]) + 100., 200.)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 2], 'ko-')
ax = p.text(0.95, 100.0, 'D', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('Number of clonal strains', fontsize = AxisLabelFontSize)
p.xlabel('Turbulence level ($ T $)', fontsize = AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks_04, size=AxisTickFontSize)
p.grid(True)

p.show()

## Printing stuff just in case...
#mutations = p.genfromtxt('MutationRecords.dat')
#MUTT = MutationPerClone(mutations)
#print "Number of cells :", MUTT[0]
#print "Number of clonal strains :", MUTT[1]
#print "---------------------------------------------"
#print "Average number of point mutations:", MUTT[2], "+/-", MUTT[3], \
#"(mean +/- STD)"
#print "Average number of duplications   :", MUTT[4], "+/-", MUTT[5], \
#"(mean +/- STD)"
#print "Average number of deletions      :", MUTT[6], "+/-", MUTT[7], \
#"(mean +/- STD)"
