#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:12:02 2023

@author: 'Denise Lavezzari'
"""

## The script takes as input the results for c-sgRNAs and creates the plot on the two replicates.

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import import_data as imp_d


# Open folder with software results taking into account the two time point replicates
os.chdir(os.getcwd()) 
rep1 = sys.argv[1]
rep2 = sys.argv[2]

name = rep1.split("_1",1)[0] + '_' + rep1.split("_",3)[2]

os.chdir(rep1)
# periscope
rep1_periscope = imp_d.periscope_import(rep1)
# leTRS
rep1_leTRS = imp_d.leTRS_import(rep1)
# sgDI-tector
rep1_sgDI = imp_d.sgDI_import(rep1)
# calculate percentage w.r.t. total fragments
rep1_mergione = imp_d.mergione(rep1_periscope, rep1_leTRS, rep1_sgDI)

os.chdir("../")
os.chdir(rep2)
# periscope
rep2_periscope = imp_d.periscope_import(rep2)
# leTRS
rep2_leTRS = imp_d.leTRS_import(rep2)
# sgDI-tector
rep2_sgDI = imp_d.sgDI_import(rep2)
# calculate percentage w.r.t. total fragments
rep2_mergione = imp_d.mergione(rep2_periscope, rep2_leTRS, rep2_sgDI)

# #Calculate mean between the two replicates
os.chdir("../")
mergione = pd.concat((rep1_mergione, rep2_mergione))
mergione = mergione.groupby('orf', sort=False).mean()
mergione = mergione.reset_index()

# create bar plot based on the mean
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
ax.set_title(name)
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector'), loc='upper left')
# plt.show()

plt.savefig(name + '_barplot_canonical.jpeg')
plt.close()

# save means to use for orf plot
mergione.to_csv(name +'_mean_orfs.csv', index = False, header=True, sep = "\t")
