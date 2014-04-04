#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates and plots the correlation between the numbers of different types of
mutations and HGT events using orthogonal distance regression (ODR). Calculates
and plots also correlation between the number of all mutations and the number of
the clonal strains using ODR.
This is HGT version (modified on 04/04/2012)

Created on Tue Jan 31 17:30:17 2012
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com

"""
import pylab as p
import scipy.odr as oder
import scipy.stats as stat
import Y_axis_div as div # this one is mine

AxisLabelFontSize = 22
AxisTickFontSize = 22
AnnotateFontSize = 19

SectOnYaxis = 7
SectOnYaxis_2 = 7
SectOnYaxis_3 = 5

#--- importing data file ---
TheFinalArray = p.genfromtxt("MutationStats.dat")
#--- defining funtions for linear regression ---
fitfunc = lambda p, x: p[0] + p[1] * x
the_line = lambda x, a, b: a * x + b
mymodel = oder.Model(fitfunc)

sum_of_mut = TheFinalArray[:, 3] + TheFinalArray[:, 5] + TheFinalArray[:, 7]
std_of_total_mut = p.sqrt(TheFinalArray[:, 4]**2 + TheFinalArray[:, 6]**2 \
    + TheFinalArray[:, 8]**2)
max_modif = (TheFinalArray[:, 3] + TheFinalArray[:, 4]).max()
max_dupl = (TheFinalArray[:, 5] + TheFinalArray[:, 6]).max()
max_del = (TheFinalArray[:, 7] + TheFinalArray[:, 8]).max()
max_hgt = (TheFinalArray[:, 9] + TheFinalArray[:, 10]).max()
max_of_maxes = max(max_modif, max_dupl, max_del, max_hgt)
[MIN, MaxMutationNumber, interv] = div.divide_Y_axis(0.0, max_of_maxes,
    SectOnYaxis)
tick_marks = p.arange(MIN, MaxMutationNumber + 0.001, interv)
Letter_loc = MaxMutationNumber - interv / 2.0
p.figure(1, figsize=(1000, 600))
#--- fitting the regression line ---
mydata_1 = oder.RealData(TheFinalArray[:, 3], TheFinalArray[:, 5],
                         sx=TheFinalArray[:, 4], sy=TheFinalArray[:, 6])
myodr_1 = oder.ODR(mydata_1, mymodel, beta0=[1., 1.], maxit=1000)
myoutput_1 = myodr_1.run()
print "# Duplications vs modifications"
myoutput_1.pprint()
per_1 = stat.spearmanr(TheFinalArray[:, 3], TheFinalArray[:, 5])
#--- and plotting... ---
p.subplot(221)
p.plot(TheFinalArray[:, 3], the_line(TheFinalArray[:, 3], myoutput_1.beta[1],
       myoutput_1.beta[0]), 'k-')
p.errorbar(TheFinalArray[:, 3], TheFinalArray[:, 5], TheFinalArray[:, 6],
           TheFinalArray[:, 4], 'ko')
ax = p.text(Letter_loc, Letter_loc, 'A', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
ax = p.text(MIN + interv / 2.0, Letter_loc, r'$\rho = %1.2f$'%(per_1[0],),
            horizontalalignment='center', verticalalignment='center',
            fontweight='bold', fontsize=23)
p.axis([0.0, MaxMutationNumber, 0.0, MaxMutationNumber])
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.xlabel('number of modifications', fontsize = AxisLabelFontSize)
p.ylabel('number of duplications', fontsize = AxisLabelFontSize)
p.grid(True)

#--- fitting the regression line ---
mydata_2 = oder.RealData(TheFinalArray[:, 3], TheFinalArray[:, 7],
                         sx=TheFinalArray[:, 4], sy=TheFinalArray[:, 8])
myodr_2 = oder.ODR(mydata_2, mymodel, beta0=[1., 1.], maxit=1000)
myoutput_2 = myodr_2.run()
print "# Deletions vs modifications"
myoutput_2.pprint()
per_2 = stat.spearmanr(TheFinalArray[:, 3], TheFinalArray[:, 7])
#--- and plotting... ---
p.subplot(222)
p.plot(TheFinalArray[:, 3], the_line(TheFinalArray[:, 3], myoutput_2.beta[1],
       myoutput_2.beta[0]), 'k-')
p.errorbar(TheFinalArray[:, 3], TheFinalArray[:, 7], TheFinalArray[:, 8],
            TheFinalArray[:, 4], 'ko')
ax = p.text(Letter_loc, Letter_loc, 'B', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
ax = p.text(MIN + interv / 2.0, Letter_loc, r'$\rho = %1.2f$'%(per_2[0],),
            horizontalalignment='center', verticalalignment='center',
            fontweight='bold', fontsize=23)
p.axis([0.0, MaxMutationNumber, 0.0, MaxMutationNumber])
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.xlabel('number of modifications', fontsize = AxisLabelFontSize)
p.ylabel('number of deletions', fontsize = AxisLabelFontSize)
p.grid(True)

#--- fitting the regression line ---
mydata_3 = oder.RealData(TheFinalArray[:, 5], TheFinalArray[:, 7],
                         sx=TheFinalArray[:, 6], sy=TheFinalArray[:, 8])
myodr_3 = oder.ODR(mydata_3, mymodel, beta0=[1., 1.], maxit=1000)
myoutput_3 = myodr_3.run()
print "# Deletions vs duplications"
myoutput_3.pprint()
per_3 = stat.spearmanr(TheFinalArray[:, 5], TheFinalArray[:, 7])
#--- and plotting... ---
p.subplot(223)
p.plot(TheFinalArray[:, 5], the_line(TheFinalArray[:, 5], myoutput_3.beta[1],
       myoutput_3.beta[0]), 'k-')
p.errorbar(TheFinalArray[:, 5], TheFinalArray[:, 7], TheFinalArray[:, 8],
           TheFinalArray[:, 6], 'ko')
ax = p.text(Letter_loc, Letter_loc, 'C', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
ax = p.text(MIN + interv / 2.0, Letter_loc, r'$\rho = %1.2f$'%(per_3[0],),
            horizontalalignment='center', verticalalignment='center',
            fontweight='bold', fontsize=23)
p.axis([0.0, MaxMutationNumber, 0.0, MaxMutationNumber])
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.xlabel('number of duplicatons', fontsize = AxisLabelFontSize)
p.ylabel('number of deletions', fontsize = AxisLabelFontSize)
p.grid(True)

#--- fitting the regression line ---
mydata_5 = oder.RealData(TheFinalArray[:, 5], TheFinalArray[:, 9],
                         sx=TheFinalArray[:, 6], sy=TheFinalArray[:, 10])
myodr_5 = oder.ODR(mydata_5, mymodel, beta0=[1., 1.], maxit=1000)
myoutput_5 = myodr_5.run()
print "# Deletions vs HGT events"
myoutput_5.pprint()
per_5 = stat.spearmanr(TheFinalArray[:, 5], TheFinalArray[:, 9])
#--- and plotting... ---
p.subplot(224)
p.plot(TheFinalArray[:, 5], the_line(TheFinalArray[:, 5], myoutput_5.beta[1],
       myoutput_5.beta[0]), 'k-')
p.errorbar(TheFinalArray[:, 5], TheFinalArray[:, 9], TheFinalArray[:, 10],
           TheFinalArray[:, 6], 'ko')
ax = p.text(Letter_loc, Letter_loc, 'D', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
ax = p.text(MIN + interv / 2.0, Letter_loc, r'$\rho = %1.2f$'%(per_5[0],),
            horizontalalignment='center', verticalalignment='center',
            fontweight='bold', fontsize=23)
p.axis([0.0, MaxMutationNumber, 0.0, MaxMutationNumber])
p.xticks(size=AxisTickFontSize)
p.yticks(tick_marks, size=AxisTickFontSize)
p.xlabel('number of duplicatons', fontsize = AxisLabelFontSize)
p.ylabel('number of HGT events', fontsize = AxisLabelFontSize)
p.grid(True)


#--- fitting the regression line ---
mydata_4 = oder.RealData(TheFinalArray[:, 2], sum_of_mut, sy=std_of_total_mut)
myodr_4 = oder.ODR(mydata_4, mymodel, beta0=[1., 1.], maxit=1000)
myoutput_4 = myodr_4.run()
print "# All mutations vs clone number:"
myoutput_4.pprint()
per_4 = stat.spearmanr(TheFinalArray[:, 2], sum_of_mut)
#--- and plotting... ---
p.figure(2,  figsize=(1000, 600))
max_clone_num = max(TheFinalArray[:, 2])
[MIN, MaxClonalStrain, intervvv] = div.divide_Y_axis(0.0, max_clone_num,
    SectOnYaxis_2)
tick_marks_2 = p.arange(MIN, MaxClonalStrain + 0.1, intervvv)
[MIN, SummOfMut, intervvvvv] = div.divide_Y_axis(0.0, sum_of_mut.max(),
    SectOnYaxis_3)
tick_marks_3 = p.arange(MIN, SummOfMut + 0.1, intervvvvv)
Letter_loc_2 = MaxClonalStrain - intervvv / 2.0
Letter_loc_3 = SummOfMut - intervvvvv / 2.0
p.plot(TheFinalArray[:, 2], the_line(TheFinalArray[:, 2], myoutput_4.beta[1],
       myoutput_4.beta[0]), 'k-')
p.errorbar(TheFinalArray[:, 2], sum_of_mut, std_of_total_mut, fmt='ko',
           ecolor='k')
ax = p.text(Letter_loc_2, Letter_loc_3, 'E', horizontalalignment='center',
            verticalalignment='center', fontweight='bold', fontsize=23)
ax = p.text(MIN + intervvv, Letter_loc_3, r'$\rho = %1.2f$'%(per_4[0],),
            horizontalalignment='center', verticalalignment='center',
            fontweight='bold', fontsize=23)
p.xlim(xmin=0)
p.ylim(ymin=0)
p.xticks(tick_marks_2, size=AxisTickFontSize)
p.yticks(tick_marks_3, size=AxisTickFontSize)
p.xlabel('number of clonal strains', fontsize = AxisLabelFontSize)
p.ylabel('number of all mutations', fontsize = AxisLabelFontSize)
p.grid(True)

p.show()