#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:17:52 2022

@author: 'Denise Lavezzari'
"""

# The script for each sample takes as input the results for c-sgRNAs obtained by the three software and compare it creating Venn and barplots.

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted
import import_data as imp_d

# Open folder with software results
os.chdir(os.getcwd()) 
sample = sys.argv[1]

os.chdir(sample)

# orf dataframe
orfs = pd.DataFrame()
orfs['orf'] = ['ORF 1ab','ORF S','ORF 3a','ORF E','ORF M','ORF 6', 'ORF 7a', 'ORF 7b', 'ORF 8', 'ORF 10']

# periscope
sample_periscope = imp_d.periscope_import(sample)
sub_periscope = pd.DataFrame()
sub_periscope[['orf', 'counts']] = sample_periscope[['orf', 'sgRNA_count']]
df1 = set(sub_periscope['orf'].tolist())

# leTRS
sample_leTRS = imp_d.leTRS_import(sample)
sub_leTRS = pd.DataFrame()
sub_leTRS[['orf', 'counts']] = sample_leTRS[['subgenome', 'peak_count']]
df2 = set(sub_leTRS['orf'].tolist())

# sgDI-tector
sample_sgDI = imp_d.sgDI_import(sample)
sub_sgDI = pd.DataFrame()
sub_sgDI[['orf', 'counts']] = sample_sgDI[['ORF', 'counts']]
df3 = set(sub_sgDI['orf'].tolist())

# Create Venn diagram to compare the results 
plt.figure(figsize=(4,4), dpi=600)
venn3_unweighted([df1, df2, df3], set_labels = ('Periscope','LeTRS','sgDI-Tector'), set_colors=('darkgreen','blue','red'), alpha = 0.4)
plt.savefig(sample + '_venn3_canonical.jpeg')
plt.close()

# calculate concordance rate
cr = round(len(set.intersection(df1, df2, df3)) / len(set.union(df1, df2, df3)) * 100, 2)
data = pd.DataFrame({'sample': sample, 'concordance_rate': cr}, index=[0])
data.to_csv(r'concordance_percentage_canonical.csv', index = False, header=True, sep = "\t")

# calculate percentage w.r.t. total fragments
mergione = orfs.merge(sub_periscope, on='orf', how='outer').merge(sub_leTRS,on='orf', how='outer', suffixes=('_periscope', '_leTRS')).merge(sub_sgDI,on='orf', how='outer').fillna(0)

mergione['%periscope'] = round((mergione['counts_periscope'] / 421872) * 100,5)
mergione['%leTRS'] = round((mergione['counts_leTRS'] / 421872) * 100,5)
mergione['%sgDI-Tector'] = round((mergione['counts'] / 421872) * 100,5)

mergione.to_csv(r'fragments_counts.csv', index = False, header=True, sep = "\t")

# create bar plot
fig, (ax, ax2, ax3) = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [4, 3, 3]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(mergione))
width = 0.2

bar1 = ax.bar(x_axis, mergione['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2= ax.bar(x_axis + width, mergione['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, mergione['counts'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, mergione['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, mergione['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, mergione['counts'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar7 = ax3.bar(x_axis, mergione['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar8 = ax3.bar(x_axis + width, mergione['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar9 = ax3.bar(x_axis + width*2, mergione['counts'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax3.set_ylim(0, 10)
ax2.set_ylim(10, 100)
ax.set_ylim(100,8000)  # outliers only

# hide the spines between ax and ax2
ax3.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)
ax3.tick_params(labeltop=False, bottom=False)

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax3.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax3.set_xticks(x_axis+width,mergione['orf'], fontsize=7, rotation=90)
ax2.set_ylabel('Counts', loc='center')
ax.set_title(sample.split("_NT",1)[0])
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector'), loc='upper left')
plt.savefig(sample + '_barplot_canonical.jpeg')
plt.close()


