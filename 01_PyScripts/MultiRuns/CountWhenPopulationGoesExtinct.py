#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 16:17:53 2010

@author: piotr bentkowski
p.bentkowski@uea.ac.uk
"""
#import pylab as p
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import linecache as ln


#--- read data from files in directiories ---
def LoadEnvelopeMeanSize(arg, dirname, files):
    for file in files:
        filepath = os.path.join(dirname, file)
        if filepath == os.path.join(dirname, 'GeneralData.dat'):
            data = np.genfromtxt(filepath)
            counter = np.nan
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


def main():
    """ """
    TheData = []
    os.path.walk(os.getcwd(), LoadEnvelopeMeanSize, TheData)
    TheData = sorted(TheData, key=lambda data_sort: data_sort[1])

    #  ---- loading and cheking data ---
    Die = np.zeros((len(TheData), 2))
    j = 0
    for item in TheData:
        Die[j, 0] = item[0]
        Die[j, 1] = item[1]
        j += 1

    Die = Die[~np.isnan(Die).any(1)]  # removing NaN's from the output
#    print Die
    np.savetxt("Extinct.dat", Die, fmt='%1.5f')

    print len(Die), Die[0][1]
    # === trim die ===
    tmpDie = []
    check = Die[0][1]
    for ii, itm in enumerate(Die):
        oneDie = []
        for jj in np.arange(ii+1, len(Die)):
            if itm[1] == Die[jj][1]:
                oneDie.append(Die[jj][0])
        if itm[1] != check:
            oD = np.array(oneDie)
            tmpDie.append((itm[1], np.mean(oD), np.std(oD)))
            check = itm[1]
    Die = np.array(tmpDie)

    print ""
    print Die

    AxisLabelFontSize = 18
    AxisTickFontSize = 17

    plt.figure(9, figsize=(10, 6))
    plt.plot(Die[:, 0], Die[:, 1], 'ko-')
#    plt.yscale('log', basey=10)
    plt.fill_between(Die[:, 0], Die[:, 1] + Die[:, 2], Die[:, 1] - Die[:, 2],
                     color=(0.75, 0.75, 0.75, 0.75))
    plt.xlabel('random death rate ($ \delta $)', fontsize=AxisLabelFontSize)
    plt.ylabel('time when population goes extinct (steps)',
               fontsize=AxisLabelFontSize)
    plt.axis([0, 0.25, 0, 25])
    plt.xticks(size=AxisTickFontSize)
    plt.yticks(size=AxisTickFontSize)
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("figure_S1.pdf", dpi=1000, pad_inches=0.0)
    plt.show()

if __name__ == "__main__":
    main()
