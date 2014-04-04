#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reads a few files and calculates what is the relation between available resources
and the gene number that the population of cells can maintain.

Genome Streamlining Project.
Created on Thu Jan 20 18:20:15 2011
Author: Piotr Bentkowski - p.bentkowski@uea.ac.uk, bentkowski.piotr@gmail.com
"""

import pylab as p
import re
import linecache as ln

FontSize = 17
TickSize = 14
#LineSize = 10

data = p.genfromtxt("GeneralData.dat")
if (data[-1,2] == 0.0):
    print "This ecosystem died out before the end of the simulation. Sorry :-("
    exit()
par_0 = re.split(" ", ln.getline('ModelParams.dat', 6))
tubulence = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 21))
core_genes = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 22))
gene_cost = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 23))
max_uptake_value = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 24))
const_cost = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 25))
minimal_size = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 26))
reproduction_size = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 29))
death_rate = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 45))
ASLD_hist_ranges_size = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 46))
ASLD_hist_ranges_bin_width = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 49))
uptake_repr_hist = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 50))
uptake_repr_bin = float(par_0[6])
par_0 = re.split(" ", ln.getline('ModelParams.dat', 51))
gene_num_bins = float(par_0[6])

ASLD = p.genfromtxt("AgeTillReproduction.dat")
tot_gain_at_repr = p.genfromtxt("AvarGainAtRepr.dat")
intake_at_repr = p.genfromtxt("AvarIntakeAtRepr.dat")
genome_size = p.genfromtxt("GeneNumberOfRepr.dat")

# cropping data to the second, most interesting half
ASLD = ASLD[-p.floor(ASLD.shape[0]/2):]
tot_gain_at_repr = tot_gain_at_repr[-p.floor(tot_gain_at_repr.shape[0]/2):]
intake_at_repr = intake_at_repr[-p.floor(intake_at_repr.shape[0]/2):]
data = data[-p.floor(data.shape[0]/2):]
genome_size = genome_size[-p.floor(genome_size.shape[0]/2):]

N_max = 60.0
x = p.arange(0.0, N_max + 1.0, 1.0)
y = gene_cost * (x + x**2 + core_genes) + const_cost
Y_max = p.ceil(y.max() / 10.0) * 10.0
real_mean_number_of_genes = data[:,2].mean()
real_gene_cost = gene_cost * (real_mean_number_of_genes * \
    (1.0 + real_mean_number_of_genes) + core_genes) + const_cost

uptake_hist_ranges = p.arange(0.0, uptake_repr_hist * uptake_repr_bin,
                             uptake_repr_bin)
ASLD_hist_ranges = p.arange(0.0, ASLD_hist_ranges_size \
                * ASLD_hist_ranges_bin_width, ASLD_hist_ranges_bin_width)
genome_size_bins = p.arange(1.0, gene_num_bins + 1, 1.0)

#print tot_gain_at_repr.shape[0], tot_gain_at_repr.shape[1]
#print uptake_hist_ranges.shape[0], ASLD_hist_ranges.shape[0]
#print ASLD[1, :].shape[0]

tot_gain_means = p.zeros(tot_gain_at_repr.shape[0])
intake_means = p.zeros(tot_gain_at_repr.shape[0])
genome_size_means = p.zeros(tot_gain_at_repr.shape[0])
ASLD_means = p.zeros(tot_gain_at_repr.shape[0])

for i in xrange(0, tot_gain_at_repr.shape[0], 1):
    ASLD_means[i] = p.sum(ASLD[i, :] * ASLD_hist_ranges ) / p.sum(ASLD[i, :])
#    ASLD[i, :] = ASLD[i, :] / ASLD[i, :].sum()
    tot_gain_means[i] = p.sum(tot_gain_at_repr[i, :] * uptake_hist_ranges) \
        / p.sum(tot_gain_at_repr[i, :])
    intake_means[i] = p.sum(intake_at_repr[i, :] * uptake_hist_ranges) \
        / p.sum(intake_at_repr[i, :])
    genome_size_means[i] = p.sum(genome_size[i, :] * genome_size_bins) \
        /p.sum(genome_size[i, :])

Summ_ASLD = p.zeros(ASLD.shape[1])
ASLD_to_hist_mean = p.zeros(ASLD.shape[1])
ASLD_to_hist_std = p.zeros(ASLD.shape[1])
for j in xrange(0, ASLD.shape[1], 1):
    ASLD_to_hist_mean[j] = ASLD[:, j].mean()
    ASLD_to_hist_std[j] = ASLD[:, j].std()
    Summ_ASLD[j] = p.sum(ASLD[:, j])

#Summ_ASLD = Summ_ASLD[~p.isnan(Summ_ASLD).any(0)]
#The_summ_ASLD = Summ_ASLD.mean()
#ASLD_to_hist_mean = ASLD_to_hist_mean / Summ_ASLD
#ASLD_to_hist_std = ASLD_to_hist_std / Summ_ASLD

tot_gain_means = tot_gain_means[~p.isnan(tot_gain_means).any(0)]
intake_means = intake_means[~p.isnan(intake_means).any(0)]
genome_size_means =  genome_size_means[~p.isnan(genome_size_means).any(0)]

repr_mean_number_of_genes = genome_size_means.mean()
real_gene_cost = gene_cost * (repr_mean_number_of_genes \
    + repr_mean_number_of_genes**2 + core_genes) + const_cost

MEAN_tot_gain = tot_gain_means.mean()
MEAN_itake = intake_means.mean()
The_diff = tot_gain_means - intake_means
MEAN_diff = The_diff.mean()
ASLD_means = ASLD_means[~p.isnan(ASLD_means).any(0)]
ASLD_tot_mean = ASLD_means.mean()

print "ASLD = the age since the last division"
print "-----------------------------------------------------------"
print "The mean cells total gain            :", MEAN_tot_gain
print "The mean cells net uptake            :", MEAN_itake
print "The difference of gain - uptake      :", MEAN_diff
print "The real cost of living of repr cells:", real_gene_cost
print "-----------------------------------------------------------"
print "The mean ASLD          :", ASLD_tot_mean
print "Minimal permitted ASLD :", p.floor((reproduction_size * 0.5) \
/ max_uptake_value)
print "Life expectancy        :", 1.0 / (data[:, 11] / data[:, 4]).mean()
print "-----------------------------------------------------------"
print "Maximal permitted net uptake"
print "  with this mean number of genes: ", max_uptake_value \
- real_gene_cost
print "-----------------------------------------------------------"
print "Cheking (mean ASLD) x (net uptake):", ASLD_tot_mean * MEAN_itake

## plotting the results
p.figure(1, figsize=(1000,600))
LineWidth = 4.5
p.title('Turbulence level = %(turb)1.3f ; Death rate = %(death)1.3f'% \
        {"turb":tubulence, "death":death_rate})
p.bar(ASLD_hist_ranges, ASLD_to_hist_mean + ASLD_to_hist_std,
      width=ASLD_hist_ranges_bin_width, color=(0.75,0.75,0.75,0.75), linewidth=0)
p.bar(ASLD_hist_ranges, ASLD_to_hist_mean, width=ASLD_hist_ranges_bin_width,
      color='k', linewidth=0)
p.axis([0, ASLD_hist_ranges.max(), 0, (ASLD_to_hist_mean+ASLD_to_hist_std).max()])
p.xlabel('Time passed between reproductions (steps)', fontsize = FontSize)
#p.ylabel('Frequency of occurrence', fontsize = FontSize)
p.ylabel('Number of cells', fontsize = FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)

p.figure(2, figsize=(1000,600))
p.title('Turbulence level = %(turb)1.3f ; Death rate = %(death)1.3f'% \
        {"turb":tubulence, "death":death_rate})
p.plot(x, y, 'ko-', markersize=4)
p.title('Gene number of cells reproducting (red vertical) and all cells (blue vertical)')
p.hlines(MEAN_diff, 0.0, N_max, 'g', linestyles='dashed', lw=2)
#p.hlines(real_gene_cost, 0.0, N_max, 'k', linestyles='solid')
p.vlines(repr_mean_number_of_genes, 0.0, Y_max , 'r', linestyles='dashed', lw=2)
p.vlines(real_mean_number_of_genes, 0.0, Y_max , 'b', linestyles='dashed', lw=2)
#p.axis([0.0, N_max, 0.0, Y_max])
p.axis([0.0, 50.0001, 0.0, 8.0001])
p.xlabel('Number of metabolic genes', fontsize = FontSize)
p.ylabel('Resources', fontsize = FontSize)
p.xticks(size=TickSize)
p.yticks(size=TickSize)
p.grid(True)

p.show()
