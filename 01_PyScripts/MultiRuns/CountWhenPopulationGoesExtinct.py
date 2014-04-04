#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 16:17:53 2010

@author: piotr bentkowski
p.bentkowski@uea.ac.uk
"""
import pylab as p
import os
import re
import linecache as ln

AxisLabelFontSize = 16
AxisTickFontSize = 15


#--- read data from files in directiories ---
def LoadEnvelopeMeanSize(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GeneralData.dat'):
            data = p.genfromtxt(filepath)
            counter = p.nan
            if data[-1, 4] == 0.0:
                counter = 0.0
                for i in data[:, 4]:
                    if i != 0.0:
                        counter += 1
                    else:
                        break
                filepath_death = os.path.join(dirname, 'ModelParams.dat')
                l = re.split(" ", ln.getline(filepath_death, 29))
                death_rate = float(l[6])
                arg.append((counter, death_rate))

TheData = []
os.path.walk(os.getcwd(), LoadEnvelopeMeanSize, TheData)
TheData = sorted(TheData, key=lambda data_sort: data_sort[1])

#  ---- loading and cheking data ---
Die = p.zeros((len(TheData), 2))
j = 0
for item in TheData:
    Die[j, 0] = item[0]
    Die[j, 1] = item[1]
    j += 1

Die = Die[~p.isnan(Die).any(1)]  # removing NaN's from the output
print Die
p.savetxt("Extinct.dat", Die, fmt='%1.5f')
#Die[:,0] = p.log10(Die[:, 0])

p.figure(9, figsize=(1000, 600))
#p.get_current_fig_manager().window.wm_geometry("1000x600+200+110")
p.plot(Die[:, 1], Die[:, 0], 'ko-')
#p.yscale('log', basey = 10)
#p.plot(Die[:,2],Die[:,0]-Die[:,1], 'k--' )
#p.plot(Die[:,2],Die[:,0]+Die[:,1], 'k--' )
#p.errorbar(Die[:,2],Die[:,0], yerr=Die[:,1], fmt='o-', ecolor='k')
p.xlabel('death rate ($ \delta $)', fontsize=AxisLabelFontSize)
p.ylabel('time when population goes extinct (steps)',
         fontsize=AxisLabelFontSize)
#p.axis([0, Die[:, 1].max() + 0.05 * Die[:, 1].max(), 0,
#        Die[:, 0].max() + 0.02 * Die[:, 0].max()])
p.xticks(size=AxisTickFontSize)
p.yticks(size=AxisTickFontSize)
p.grid(True)
p.show()
