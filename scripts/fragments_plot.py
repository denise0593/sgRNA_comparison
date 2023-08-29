#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:40:28 2023

@author: 'Denise Lavezzari'
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

os.chdir("Subgenomic_RNA/")

####################################################################
######################### CANONICAL ################################
####################################################################

# Prepare data for the analysis: calculate the percentage, add the column specifing the software ans simplify sample name

df_per_c = pd.read_excel('Suppl. Table S2A.xls', skiprows=5,  usecols=[0,10], header=1, names=['Sample','% fragments'])
df_per_c["% fragments"] = df_per_c["% fragments"]*100
df_per_c["Software"] = "Periscope"
df_per_c['Sample'] = df_per_c['Sample'].str.split("_").str[0]  + '_' + df_per_c['Sample'].str.split("_").str[2]

df_leTRS_c = pd.read_excel('Suppl. Table S2A.xls', skiprows=5,  usecols=[0,11], header=1, names=['Sample', '% fragments']) 
df_leTRS_c["% fragments"] = df_leTRS_c["% fragments"]*100
df_leTRS_c["Software"] = "LeTRS"
df_leTRS_c['Sample'] = df_leTRS_c['Sample'].str.split("_").str[0] + '_' + df_leTRS_c['Sample'].str.split("_").str[2]

df_sgDI_c = pd.read_excel('Suppl. Table S2A.xls', skiprows=5,  usecols=[0,12], header=1, names=['Sample','% fragments']) 
df_sgDI_c["% fragments"] = df_sgDI_c["% fragments"]*100
df_sgDI_c["Software"] = "sgDI-tector"
df_sgDI_c['Sample'] = df_sgDI_c['Sample'].str.split("_").str[0] + '_' + df_sgDI_c['Sample'].str.split("_").str[2]

full_df_c = pd.concat((df_per_c,df_leTRS_c,df_sgDI_c))

### PLOT #######

plt.figure(dpi=600)
a = sns.barplot(x="Sample", y="% fragments",  hue="Software", data=full_df_c, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black')
a = sns.swarmplot(x="Sample", y="% fragments", hue="Software", data=full_df_c, palette='dark:0', dodge=True, size=4, edgecolor="gray", legend=False)
a.set_ylim(0, 3)
a.set(xlabel=None)

plt.savefig('plot_percentage_fragments_canonical.jpeg')
plt.close()

####################################################################
####################### NON CANONICAL ##############################
####################################################################

# Prepare data for the analysis: calculate the percentage, add the column specifing the software ans simplify sample name

df_per_nc = pd.read_excel('Suppl. Table S2A.xls', skiprows=5,  usecols=[0,16], header=1, names=['Sample','% fragments']) 
df_per_nc["% fragments"] = df_per_nc["% fragments"]*100
df_per_nc["Software"] = "Periscope"
df_per_nc['Sample'] = df_per_nc['Sample'].str.split("_").str[0] + '_' + df_per_nc['Sample'].str.split("_").str[2]


df_leTRS_nc = pd.read_excel('Suppl. Table S2A.xls', skiprows=5,  usecols=[0,17], header=1, names=['Sample','% fragments'])
df_leTRS_nc["% fragments"] = df_leTRS_nc["% fragments"]*100
df_leTRS_nc["Software"] = "LeTRS"
df_leTRS_nc['Sample'] = df_leTRS_nc['Sample'].str.split("_").str[0]+ '_' + df_leTRS_nc['Sample'].str.split("_").str[2]


df_sgDI_nc = pd.read_excel('Suppl. Table S2A.xls', skiprows=5,  usecols=[0,18], header=1, names=['Sample','% fragments']) 
df_sgDI_nc["% fragments"] = df_sgDI_nc["% fragments"]*100
df_sgDI_nc["Software"] = "sgDI-tector"
df_sgDI_nc['Sample'] = df_sgDI_nc['Sample'].str.split("_").str[0] + '_' + df_sgDI_nc['Sample'].str.split("_").str[2]


full_df_nc = pd.concat((df_per_nc,df_leTRS_nc,df_sgDI_nc))


### PLOT #######

fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [7, 3]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(full_df_nc['Sample'].unique()))
width = 0.2

sns.barplot(x="Sample", y="% fragments",  hue="Software", data=full_df_nc, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black', ax=axes[0])
sns.swarmplot(x="Sample", y="% fragments", hue="Software", data=full_df_nc, palette='dark:0', dodge=True, size=4, edgecolor="gray", legend=False, ax=axes[0])

sns.barplot(x="Sample", y="% fragments",  hue="Software", data=full_df_nc, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black', ax=axes[1])
sns.swarmplot(x="Sample", y="% fragments", hue="Software", data=full_df_nc, palette='dark:0', dodge=True, size=4, edgecolor="gray", legend=False, ax=axes[1])


# remove legend repetition
axes[1].get_legend().set_visible(False)

# remove y label repetition an set y label
axes[1].set(ylabel=None) 
axes[0].set(ylabel=None) 
axes[0].set_ylabel('% fragments', loc='center')

# remove x label repetition and set the one on the bottom subplot
axes[0].set(xlabel=None) 
axes[1].set(xlabel=None) 
axes[0].set_xticks([])

axes[1].set_ylim(0, 0.05)
axes[0].set_ylim(0.5,3)  # outliers only

# hide the spines between ax and ax2
axes[1].spines['top'].set_visible(False)
axes[0].spines['bottom'].set_visible(False)

# set diagonal tiks 
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False)
axes[1].plot((-d, +d), (-d, +d), **kwargs)         # bottom-left diagonal 
axes[1].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-right diagonal 

kwargs.update(transform=axes[1].transAxes)  # switch to the bottom axes
axes[0].plot((-d, +d), (1 - d, 1 + d), **kwargs) # top-left diagonal
axes[0].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # top-right diagonal

plt.savefig('plot_percentage_fragments_non_canonical.jpeg')
plt.close()
