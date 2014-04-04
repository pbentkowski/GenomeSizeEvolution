#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the avarage number of born and dead cells of the stable part of the
run comparing GeneralData.dat in different directories (between different runs).

Created on Fri Jan 13 20:10:18 2012

Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import os, re
import linecache as ln
# the one below is mine
import file_len as fl

#--- stuff which makes the figure look nice
TickFont = 16
LabelFont = 20
TextFontSize = 20
TextPos_X = 0.475
TextPos_Y = 35.0
Y_min = 0
Y_max = 40
Y_ticks_range = p.arange(Y_min, Y_max + 0.0001, 10)

#--- read data from files in directiories ---
def LoadTheBigData(arg, dirname, files):
    '''Walks thorough directory tree and looks into GeneralData.dat files
    calculating the average number and STD of born and dead cells per time step,
    but only if population makes it to the end of simulation. Appends all of
    that to a list'''
    for file in files:
        filepath_general = os.path.join(dirname, file)
        if filepath_general == os.path.join(dirname,'GeneralData.dat'):
            l = re.split(" ", ln.getline(filepath_general,
                                         fl.file_len(filepath_general)-1))
            is_alive = float(l[4])
            if  is_alive > 0.0:
                data = p.genfromtxt(filepath_general)
                filepath_turb = os.path.join(dirname,'ModelParams.dat')
                l = re.split(" ", ln.getline(filepath_turb, 6))
                turb_param = float(l[6])
                half_run_len = p.floor( data.shape[0] / 2.0)
                data = data[half_run_len:-1, :]
                # goes: turbulence, mean_born, std_born, mean_dead, std_dead
                arg.append((turb_param, data[:, 10].mean(), data[:, 10].std(),
                           data[:, 11].mean(), data[:, 11].std() ))
            else:
                break


TheBigData = []
os.path.walk(os.getcwd(), LoadTheBigData, TheBigData)
TheBigData = sorted(TheBigData, key=lambda data_sort: data_sort[0])

DATA = p.zeros((len(TheBigData), 5))
i = 0
for item in TheBigData:
    print "For turbulence", item[0], ":"
    print "born:", item[1], "+/-", item[2], "(mean +/- STD)"
    print "dead:", item[3], "+/-", item[4], "(mean +/- STD)"
    DATA[i, 0] = item[0]
    DATA[i, 1] = item[1]
    DATA[i, 2] = item[2]
    DATA[i, 3] = item[3]
    DATA[i, 4] = item[4]
    i += 1

p.figure(1, figsize=(1000,600))
p.subplot(211)
p.fill_between(DATA[:, 0], DATA[:, 1] + DATA[:, 2], DATA[:, 1] - DATA[:, 2],
               color=(0.75, 0.75, 0.75, 0.75))
p.plot(DATA[:, 0],DATA[:, 1], 'ko-')
p.text(TextPos_X, TextPos_Y, 'A', horizontalalignment='center',
       verticalalignment='center', fontweight='bold', fontsize=TextFontSize)
#p.xlabel('Turbulence level ($ T $)', fontsize=16)
p.ylabel('number of born', fontsize=LabelFont)
p.xticks(size=TickFont)
p.yticks(Y_ticks_range, size=TickFont)
p.ylim( Y_min, Y_max )
p.xlim(0, 0.5)
p.grid(True)

p.subplot(212)
p.fill_between(DATA[:, 0], DATA[:, 3] + DATA[:, 4], DATA[:, 3] - DATA[:, 4],
               color=(0.75, 0.75, 0.75, 0.75))
p.plot(DATA[:, 0], DATA[:, 3], 'ko-')
p.text(TextPos_X, TextPos_Y, 'B', horizontalalignment='center',
       verticalalignment='center', fontweight='bold', fontsize=TextFontSize)
p.xlabel('turbulence level ($ T $)', fontsize=LabelFont)
p.ylabel('number of dead', fontsize=LabelFont)
p.xticks(size=TickFont)
p.yticks(Y_ticks_range, size=TickFont)
p.ylim( Y_min, Y_max )
p.xlim(0, 0.5)
p.grid(True)

p.show()