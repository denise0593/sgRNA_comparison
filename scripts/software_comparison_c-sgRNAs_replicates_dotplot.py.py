#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:28:11 2023

@author: 'Denise Lavezzari'
"""

import os

os.chdir("scripts")

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import import_data as imp_d
import seaborn as sns

os.chdir("NT")

rep1 = "6265_1_24h_NT"
rep2 = "6265_2_24h_NT"

# os.chdir(os.getcwd()) 
# rep1 = sys.argv[1]
# rep2 = sys.argv[2]


orfs = pd.DataFrame()
orfs['orf'] = ['ORF 1ab','ORF S','ORF 3a','ORF E','ORF M','ORF 6', 'ORF 7a', 'ORF 7b', 'ORF 8', 'ORF 10']
 
###### REPLICATE 1 #####################
name = rep1.split("_1",1)[0] + '_' + rep1.split("_",3)[2]
os.chdir(rep1)
# periscope
rep1_periscope = imp_d.periscope_import(rep1)
rep1_periscope= orfs.merge(rep1_periscope, on='orf', how='outer').fillna(0)
rep1_periscope["Software"] = "Periscope"
# leTRS
rep1_leTRS = imp_d.leTRS_import(rep1)
rep1_leTRS= orfs.merge(rep1_leTRS, on='orf', how='outer').fillna(0)
rep1_leTRS["Software"] = "LeTRS"
# sgDI-tector
rep1_sgDI = imp_d.sgDI_import(rep1)
rep1_sgDI= orfs.merge(rep1_sgDI, on='orf', how='outer').fillna(0)
rep1_sgDI["Software"] = "sgDI-tector"
# calculate percentage w.r.t. total fragments
rep1_mergione = imp_d.mergione_vert(rep1_periscope, rep1_leTRS, rep1_sgDI)



###### REPLICATE 2 #####################
os.chdir("NT")
os.chdir(rep2)
# periscope
rep2_periscope = imp_d.periscope_import(rep2)
rep2_periscope= orfs.merge(rep2_periscope, on='orf', how='outer').fillna(0)
rep2_periscope["Software"] = "Periscope"
# leTRS
rep2_leTRS = imp_d.leTRS_import(rep2)
rep2_leTRS= orfs.merge(rep2_leTRS, on='orf', how='outer').fillna(0)
rep2_leTRS["Software"] = "LeTRS"
# sgDI-tector
rep2_sgDI = imp_d.sgDI_import(rep2)
rep2_sgDI= orfs.merge(rep2_sgDI, on='orf', how='outer').fillna(0)
rep2_sgDI["Software"] = "sgDI-tector"
# calculate percentage w.r.t. total fragments
rep2_mergione = imp_d.mergione_vert(rep2_periscope, rep2_leTRS, rep2_sgDI)

# #Calculate mean between the two replicates
os.chdir("NT")
mergione = pd.concat((rep1_mergione, rep2_mergione))



###### PLOT #####################

fig, axes = plt.subplots(3, 1, gridspec_kw={'height_ratios': [4, 3, 3]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(mergione['orf'].unique()))
width = 0.2

sns.barplot(x="orf", y="counts",  hue="Software", data=mergione, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black', ax=axes[0])
sns.swarmplot(x="orf", y="counts", hue="Software", data=mergione, palette='dark:0', dodge=True, size=4, edgecolor="gray", ax=axes[0])

sns.barplot(x="orf", y="counts",  hue="Software", data=mergione, palette=["darkgreen","blue","red"], alpha = 0.4, edgecolor = 'black', ax=axes[1])
sns.swarmplot(x="orf", y="counts", hue="Software", data=mergione, palette='dark:0', dodge=True, size=4, edgecolor="gray", ax=axes[1])

sns.barplot(x="orf", y="counts",  hue="Software", data=mergione, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black', ax=axes[2])
sns.swarmplot(x="orf", y="counts", hue="Software", data=mergione, palette='dark:0', dodge=True, size=4, edgecolor="gray", ax=axes[2])

# remove legend repetition
axes[2].get_legend().set_visible(False)
axes[1].get_legend().set_visible(False)
axes[0].get_legend().set_visible(False)

# remove y label repetition an set y label
axes[2].set(ylabel=None) 
axes[1].set(ylabel=None) 
axes[0].set(ylabel=None) 
axes[1].set_ylabel('Counts', loc='center')

# remove x label repetition and set the one on the bottom subplot
axes[0].set(xlabel=None) 
axes[1].set(xlabel=None) 
axes[2].set(xlabel=None)

axes[0].set_xticks([])
axes[1].set_xticks([])
axes[2].set_xticklabels(mergione['orf'].unique(), rotation = 90, fontsize=7)


axes[2].set_ylim(0, 10)
axes[1].set_ylim(10, 100)
axes[0].set_ylim(100,8000)  # outliers only

# hide the spines between ax and ax2
axes[2].spines['top'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['bottom'].set_visible(False)
axes[0].spines['bottom'].set_visible(False)

# set diagonal tiks 
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False)
axes[2].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
axes[2].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=axes[1].transAxes)  # switch to the bottom axes
axes[1].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
axes[1].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
axes[1].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
axes[1].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

kwargs.update(transform=axes[2].transAxes)  # switch to the bottom axes
axes[0].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
axes[0].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

# set legend
handles, labels = axes[2].get_legend_handles_labels()
del handles[0:3]
del labels[0:3]
axes[0].legend(handles, labels, loc='upper left')

fig.suptitle(name)

# plt.show()

plt.savefig(name + '_barplot_and_point_canonical.jpeg')
plt.close()
