#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Iterates over a bunch o directories and opens the files MutationRecords.dat.
Calculates the mean and the STD of numbers of mutations of all three kinds
separately. Takes under account that there are different clones and the number
of mutations is averaged over a number of clones. Writes the results to a file
MutationStats.dat.
This is HGT version (modified on 03/04/2012)

Created on Tue Jul 12 14:27:53 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""
import pylab as p
import os
import re
import linecache as ln
import file_len as fl  # this one is mine
import Y_axis_div as div  # this one is mine

AxisLabelFontSize = 20
AxisTickFontSize = 20
AnnotateFontSize = 17

SectOnYaxis = 6
SectOnYaxis_2 = 5

time_frame = 1e5


def MutationPerClone(x):
    '''Calculates average number and STD of mutations from MutationRecords.dat
    file per clone. Mut_record[0] - number of cells, [1] - number of clones,
    [2] - mean number of point mutations, [3] - STD of point mutations,
    [4] - mean number of duplications, [5] - STD  of duplications,
    [6] - mean number of deletion, [7] - STD of deletion.'''
    Mut_record = p.zeros(10)
    Mut_record[0] = x.shape[0]
    ZERRO = p.zeros((x.shape[0], 1))
    x = p.concatenate((x, ZERRO), axis=1)
    for i in xrange(0, x.shape[0]):
        if (x[i, 5] == 0):
            for j in xrange(i+1, x.shape[0]):
                if (x[j, 0] == x[i, 0]):
                    x[j, 5] = p.nan
    x = x[~p.isnan(x).any(1)]
    Mut_record[1] = x.shape[0]
    Mut_record[2] = x[:, 1].mean()
    Mut_record[3] = x[:, 1].std()
    Mut_record[4] = x[:, 2].mean()
    Mut_record[5] = x[:, 2].std()
    Mut_record[6] = x[:, 3].mean()
    Mut_record[7] = x[:, 3].std()
    Mut_record[8] = x[:, 4].mean()
    Mut_record[9] = x[:, 4].std()
    return Mut_record


def LoadMutationData(arg, dirname, files):
    '''Walks thorough directory tree and looks into MutationRecords.dat files
    calculating the average number and STD of mutations per clone, but only if
    population makes it to the end of simulation. Append all thet that to a
    list'''
    for file in files:
        filepath_general = os.path.join(dirname, file)
        if filepath_general == os.path.join(dirname, 'GeneralData.dat'):
            l = re.split(" ", ln.getline(filepath_general,
                                         fl.file_len(filepath_general)-1))
            is_alive = float(l[4])
            if is_alive > 0.0:
                filepath_turb = os.path.join(dirname, 'ModelParams.dat')
                l = re.split(" ", ln.getline(filepath_turb, 6))
                turb_param = float(l[6])
                # change the line number in ln.getline below to extract the
                # parameter you're interested in
                l = re.split(" ", ln.getline(filepath_turb, 12))
                interesting_param = float(l[6])
                filepath_mut = os.path.join(dirname, 'MutationRecords.dat')
                Mut_raw_data = p.genfromtxt(filepath_mut)
                l_2 = re.split(" ", ln.getline(filepath_mut,
                                               fl.file_len(filepath_mut)-1))
                time_stamp = float(l_2[1])
                l_3 = re.split(" ", ln.getline(filepath_mut,
                                               fl.file_len(filepath_mut)-2))
                start_point = float(l_3[1])
                arg.append((turb_param, MutationPerClone(Mut_raw_data),
                            interesting_param, time_stamp, start_point))
            else:
                break

MutationData = []
os.path.walk(os.getcwd(), LoadMutationData, MutationData)
MutationData = sorted(MutationData, key=lambda data_sort: data_sort[0])

TheFinalArray = p.zeros((len(MutationData), 15))
i = 0
for item in MutationData:
    scale_fac = (item[3] - item[4]) / time_frame
    TheFinalArray[i, 0] = item[0]    # Turbulence parameter
    TheFinalArray[i, 1] = item[1][0]  # number of cells
    TheFinalArray[i, 2] = item[1][1]  # number of clonal strains
    TheFinalArray[i, 3] = item[1][2] / scale_fac # mean number of point mutations
    TheFinalArray[i, 4] = item[1][3] / scale_fac # STD of number of point mut.
    TheFinalArray[i, 5] = item[1][4] / scale_fac # mean number of duplications
    TheFinalArray[i, 6] = item[1][5] / scale_fac # STD of number of duplications
    TheFinalArray[i, 7] = item[1][6] / scale_fac # mean number of deletion
    TheFinalArray[i, 8] = item[1][7] / scale_fac # STD of number of deletion
    TheFinalArray[i, 9] = item[1][8] / scale_fac # mean number of HGT events
    TheFinalArray[i, 10] = item[1][9] / scale_fac # STD of number of HGT events
    TheFinalArray[i, 11] = item[2]  # interesting param (usually mutarion rate)
    TheFinalArray[i, 12] = item[3]  # after how many time steps was calculated
    TheFinalArray[i, 13] = item[4] # at what time step mutation counter was zeroed
    TheFinalArray[i, 14] = time_frame  # reference time frame
    i += 1

info_str = "#turb_param num_of_cell num_of_clones mean_point_mut "\
    + "STD_point_mut mean_dupl STD_dupl mean_del STD_del interesting_param "\
    + "calculation_time zeroed_time ref_time"
p.savetxt("MutationStats.dat", TheFinalArray, fmt='%1.5f')
doc = open('MutationStats.dat', 'a')
doc.write(info_str)

#--- plotting stuff ---
max_modif = (TheFinalArray[:, 3] + TheFinalArray[:, 4]).max()
max_dupl = (TheFinalArray[:, 5] + TheFinalArray[:, 6]).max()
max_del = (TheFinalArray[:, 7] + TheFinalArray[:, 8]).max()
max_hgt = (TheFinalArray[:, 9] + TheFinalArray[:, 10]).max()
max_of_maxes = max(max_modif, max_dupl, max_del, max_hgt)
[MIN, MaxMutationNumber, interv] = div.divide_Y_axis(0.0, max_of_maxes,
                                                     SectOnYaxis)

Letter_location = MaxMutationNumber - interv / 2.0

p.figure(1, figsize=(18, 9))
tick_marks = p.arange(0., MaxMutationNumber + 0.001, interv)
p.subplot(221)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 3], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 3] - TheFinalArray[:, 4],
               TheFinalArray[:, 3] + TheFinalArray[:, 4],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'A', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of gene modifications', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim((0.0, 0.5))
p.grid(True)

p.subplot(222)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 5], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 5] - TheFinalArray[:, 6],
               TheFinalArray[:, 5] + TheFinalArray[:, 6],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'B', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of duplications', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim((0.0, 0.5))
p.grid(True)

p.subplot(223)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 7], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 7] - TheFinalArray[:, 8],
               TheFinalArray[:, 7] + TheFinalArray[:, 8],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'C', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of deletions', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim((0.0, 0.5))
p.grid(True)

p.subplot(224)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 9], 'ko-')
p.fill_between(TheFinalArray[:, 0], TheFinalArray[:, 9] - TheFinalArray[:, 10],
               TheFinalArray[:, 9] + TheFinalArray[:, 10],
               color=(0.75, 0.75, 0.75, 0.75))
ax = p.text(0.475, Letter_location, 'D', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of HGT events', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.ylim((0.0, MaxMutationNumber))
p.xlim((0.0, 0.5))
p.grid(True)

max_clone_num = max(TheFinalArray[:, 2])
[MIN, MaxClonalStrain, intervvv] = div.divide_Y_axis(0.0, max_clone_num,
                                                     SectOnYaxis_2)

p.figure(2, figsize=(18, 9))
tick_marks_04 = p.arange(0., MaxClonalStrain + 0.1, intervvv)
p.plot(TheFinalArray[:, 0], TheFinalArray[:, 2], 'ko-')
ax = p.text(0.475, intervvv / 2.0, 'E', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
p.ylabel('number of clonal strains', fontsize=AxisLabelFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=AxisLabelFontSize)
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
