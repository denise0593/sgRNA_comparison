#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:00:46 2022

@author: 'Denise Lavezzari'
"""

## The script takes as input the results for nc-sgRNAs and creates the plot on the single sample and on the two replicates.

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted

# Open folder with software results taking into account the two time point replicates
os.chdir(os.getcwd())
sample = sys.argv[1]
sample2 = sys.argv[2]

# orf dataframe
orfs = pd.DataFrame()
df_bed = pd.DataFrame()

############### SAMPLE 1 ###########################
# periscope
sample_periscope = pd.read_csv(sample + '/' + sample + '_periscope_novel_counts.csv')
sample_periscope['pos'] = sample_periscope['orf'].str.split(
    "_", 1).str.get(1).astype(str).astype(int)
sample_periscope['sample'] = sample

# leTRS
sample_leTRS = pd.read_csv(sample + '/' +'LeTRS_output_cov0/results/novel_junction.tab', sep="\t")
sample_leTRS.drop(sample_leTRS.tail(7).index, inplace=True)
sample_leTRS.insert(0, 'sample', sample)
sample_leTRS.rename(columns={"TRS_start":"pos"}, inplace=True)
sample_leTRS['nb_count'] = sample_leTRS['nb_count'].str.split("(", 1).str.get(0).astype(str).astype(int)

# sgDI-tector
sample_sgDI = pd.read_csv(sample + '/' +'sgDI_results/non_canonical_sgDI.csv', sep=",")
sample_sgDI.insert(0, 'sample', sample)
sample_sgDI.rename(columns={"ri_pos":"pos"}, inplace=True)

############### SAMPLE 2 ###########################
# periscope
sample2_periscope = pd.read_csv(sample2 + '/' + sample2 + '_periscope_novel_counts.csv')
sample2_periscope['pos'] = sample2_periscope['orf'].str.split(
    "_", 1).str.get(1).astype(str).astype(int)
sample2_periscope['sample'] = sample2

# leTRS
sample2_leTRS = pd.read_csv(sample2 + '/' +'LeTRS_output_cov0/results/novel_junction.tab', sep="\t")
sample2_leTRS.drop(sample2_leTRS.tail(7).index, inplace=True)
sample2_leTRS.insert(0, 'sample', sample2)
sample2_leTRS.rename(columns={"TRS_start":"pos"}, inplace=True)
sample2_leTRS['nb_count'] = sample2_leTRS['nb_count'].str.split("(", 1).str.get(0).astype(str).astype(int)

# sgDI-tector
sample2_sgDI = pd.read_csv(sample2 + '/' +'sgDI_results/non_canonical_sgDI.csv', sep=",")
sample2_sgDI.insert(0, 'sample', sample2)
sample2_sgDI.rename(columns={"ri_pos":"pos"}, inplace=True)

###################
orfs = pd.concat([sample_periscope[['sample', 'pos']], sample_leTRS[['sample', 'pos']], sample_sgDI[['sample', 'pos']], sample2_periscope[['sample', 'pos']], sample2_leTRS[['sample', 'pos']], sample2_sgDI[['sample', 'pos']]]).sort_values('pos')

df_bed['orf_start'] = (orfs['pos']-10).astype(int)
df_bed['orf_end'] = (orfs['pos']+10).astype(int)
df_bed.insert(0, 'chromosome', 'MN908947.3')

# create a bed file with the positions of interest
df_bed.to_csv(r'softwares.bed', sep='\t', header=False, index=False)

os.system("bedtools merge -i softwares.bed > softwares_merged_"+ sample + "_" + sample2 + ".bed")
          
# # intersect the ORFs bed
os.system("bedtools intersect -a softwares_merged_" + sample + "_" + sample2 + ".bed -b //opt/references/SARS-CoV2/ORFs.bed -wo | cut -f1,2,3,7 > intersection_ORFS_" + sample + "_" + sample2 + ".bed")

######################################################################################
######################## SAMPLE 1 CONCORDANCE RATE ###################################
######################################################################################
df_amplicon = pd.read_csv(r'softwares_merged_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end"])

df_amplicon[['periscope', 'leTRS','sgDI']] = ""
df_amplicon[['counts_periscope', 'counts_leTRS','sgDI_count']] = ""

set_periscope = set()
set_leTRS = set()
set_sgDI = set()

# count sgRNA contained in the intervals, counts of reads supporting the intervals
for i in range(len(df_amplicon)):
    start = df_amplicon.loc[i,'start']
    end = df_amplicon.loc[i,'end']
    if sample_periscope['pos'].between(start,end).any():
        df_amplicon.loc[i,'periscope'] = sample_periscope['pos'].between(start,end).sum()
        df_amplicon.loc[i,'counts_periscope'] = sample_periscope[sample_periscope['pos'].between(start,end)]['sgRNA_count'].sum()
        set_periscope.add(df_amplicon.loc[i,'start'])
    else: 
        df_amplicon.loc[i,'periscope'] = 0
        df_amplicon.loc[i,'counts_periscope'] = 0
        
    if sample_leTRS['pos'].between(start,end).any():
        df_amplicon.loc[i,'leTRS'] = sample_leTRS['pos'].between(start,end).sum()
        df_amplicon.loc[i,'counts_leTRS'] = sample_leTRS[sample_leTRS['pos'].between(start,end)]['nb_count'].sum()
        set_leTRS.add(df_amplicon.loc[i,'start'])
    else: 
        df_amplicon.loc[i,'leTRS'] = 0
        df_amplicon.loc[i,'counts_leTRS'] = 0
        
    if sample_sgDI['pos'].between(start,end).any():
        df_amplicon.loc[i,'sgDI'] = sample_sgDI['pos'].between(start,end).sum()
        df_amplicon.loc[i,'sgDI_count'] = sample_sgDI[sample_sgDI['pos'].between(start,end)]['counts'].sum()
        set_sgDI.add(df_amplicon.loc[i,'start'])
    else: 
        df_amplicon.loc[i,'sgDI'] = 0
        df_amplicon.loc[i,'sgDI_count'] = 0
        
# calculate concordance rate
df_amplicon['sum'] = df_amplicon.loc[:,['periscope','leTRS', 'sgDI']].sum(axis='columns')
df_amplicon['count'] = df_amplicon[['periscope','leTRS', 'sgDI']].gt(0).sum(axis='columns')

common = df_amplicon['count'].gt(2).sum()
total = df_amplicon['count'].gt(0).sum()
cr = round(common / total * 100, 2)
print(cr)

data = pd.DataFrame({'sample': sample, 'concordance_rate': cr}, index=[0])

######################################################################################
######################## SAMPLE 1 VENN PLOT ##########################################
######################################################################################

# create venn plot
plt.figure(figsize=(4,4), dpi=600)
venn3_unweighted([set_periscope, set_leTRS, set_sgDI], set_labels = ('Periscope','LeTRS','sgDI-Tector'), set_colors=('darkgreen','blue','red'), alpha = 0.4)
plt.title(sample.split("_NT",1)[0])
# plt.show()
plt.savefig(sample + '_venn3_non_canonical.jpeg')
plt.close()

######################################################################################
######################## SAMPLE 1 BAR PLOT ###########################################
######################################################################################

# create bar plot
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, figsize=(20,6), dpi=600)

fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(df_amplicon))
width = 0.2

bar1 = ax.bar(x_axis, df_amplicon['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2 = ax.bar(x_axis + width, df_amplicon['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, df_amplicon['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, df_amplicon['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, df_amplicon['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, df_amplicon['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax2.set_ylim(0, 10)
ax.set_ylim(10, 100)

ax2.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)

d = .005  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
ax2.set_ylabel('Counts', loc='top')

ax2.set_xticks(x_axis + width,df_amplicon['start'], fontsize=8, rotation=90)

ax.set_title(sample.split("_NT",1)[0])
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )
plt.savefig(sample + '_barplot_non_canonical.jpeg')
plt.close()


######################################################################################
######################## SAMPLE 1 CONCORDANCE RATE ORFS REGIONS ######################
######################################################################################
df_ORFs_starting = pd.read_csv(r'intersection_ORFS_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end", "ORF"])

merge = df_ORFs_starting.merge(df_amplicon, how='outer')
conte = merge[['ORF', 'counts_periscope', 'counts_leTRS','sgDI_count']].groupby('ORF').sum().reset_index()

df_ORFs_names = pd.DataFrame()
df_ORFs_names['ORF'] = ['ORF 1ab','ORF S','ORF 3a','ORF E','ORF M','ORF 6', 'ORF 7a', 'ORF 7b', 'ORF 8', 'ORF10', 'ORF N']

df_ORFs = df_ORFs_names.merge(conte, on='ORF', how='outer').fillna(0)

df_amplicon_ORFS = df_amplicon.merge(df_ORFs_starting, how='outer')
df_amplicon_ORFS.to_csv(sample + "_novel_sgRNA_softwares_table.csv", sep=';', header=True, index=False)
   
######################################################################################
######################## SAMPLE 1 BAR PLOT ORFS REGIONS #############################
######################################################################################

fig.subplots_adjust(hspace=0.05) 
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(df_ORFs))
width = 0.2

bar1 = ax.bar(x_axis, df_ORFs['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2= ax.bar(x_axis + width, df_ORFs['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, df_ORFs['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, df_ORFs['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, df_ORFs['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, df_ORFs['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax2.set_ylim(0, 10)
ax.set_ylim(10, 100)
ax2.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)


d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xticks(x_axis+width,df_ORFs['ORF'], fontsize=7, rotation=90)
ax2.set_ylabel('Counts', loc='top')
ax.set_title(sample.split("_NT",1)[0])
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )
# plt.show()
plt.savefig(sample + '_barplot_non_canonical_ORFS.jpeg')
plt.close()

######################################################################################
######################## SAMPLE 2 CONCORDANCE RATE ###################################
######################################################################################
df2_amplicon = pd.read_csv(r'softwares_merged_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end"])

df2_amplicon[['periscope', 'leTRS','sgDI']] = ""
df2_amplicon[['counts_periscope', 'counts_leTRS','sgDI_count']] = ""

set2_periscope = set()
set2_leTRS = set()
set2_sgDI = set()

# count sgRNA contained in the intervals, counts of reads supporting the intervals
for i in range(len(df2_amplicon)):
    start = df2_amplicon.loc[i,'start']
    end = df2_amplicon.loc[i,'end']
    if sample2_periscope['pos'].between(start,end).any():
        df2_amplicon.loc[i,'periscope'] = sample2_periscope['pos'].between(start,end).sum()
        df2_amplicon.loc[i,'counts_periscope'] = sample2_periscope[sample2_periscope['pos'].between(start,end)]['sgRNA_count'].sum()
        set2_periscope.add(df2_amplicon.loc[i,'start'])
    else: 
        df2_amplicon.loc[i,'periscope'] = 0
        df2_amplicon.loc[i,'counts_periscope'] = 0
        
    if sample2_leTRS['pos'].between(start,end).any():
        df2_amplicon.loc[i,'leTRS'] = sample2_leTRS['pos'].between(start,end).sum()
        df2_amplicon.loc[i,'counts_leTRS'] = sample2_leTRS[sample2_leTRS['pos'].between(start,end)]['nb_count'].sum()
        set2_leTRS.add(df2_amplicon.loc[i,'start'])
    else: 
        df2_amplicon.loc[i,'leTRS'] = 0
        df2_amplicon.loc[i,'counts_leTRS'] = 0
        
    if sample2_sgDI['pos'].between(start,end).any():
        df2_amplicon.loc[i,'sgDI'] = sample2_sgDI['pos'].between(start,end).sum()
        df2_amplicon.loc[i,'sgDI_count'] = sample2_sgDI[sample2_sgDI['pos'].between(start,end)]['counts'].sum()
        set2_sgDI.add(df2_amplicon.loc[i,'start'])
    else: 
        df2_amplicon.loc[i,'sgDI'] = 0
        df2_amplicon.loc[i,'sgDI_count'] = 0
        
# calculate concordance rate
df2_amplicon['sum'] = df2_amplicon.loc[:,['periscope','leTRS', 'sgDI']].sum(axis='columns')
df2_amplicon['count'] = df2_amplicon[['periscope','leTRS', 'sgDI']].gt(0).sum(axis='columns')

common2 = df2_amplicon['count'].gt(2).sum()
total2 = df2_amplicon['count'].gt(0).sum()
cr2 = round(common2 / total2 * 100, 2)
print(cr2)

data = pd.DataFrame({'sample': sample2, 'concordance_rate': cr2}, index=[1])
data.to_csv(r'concordance_percentage_non_canonical.csv', index = False, header=True, sep = "\t")

######################################################################################
######################## SAMPLE 2 VENN PLOT ##########################################
######################################################################################

# create venn plot
plt.figure(figsize=(4,4), dpi=600)
venn3_unweighted([set2_periscope, set2_leTRS, set2_sgDI], set_labels = ('Periscope','LeTRS','sgDI-Tector'), set_colors=('darkgreen','blue','red'), alpha = 0.4)
plt.title(sample2.split("_NT",1)[0])
# plt.show()
plt.savefig(sample2 + '_venn3_non_canonical.jpeg')
plt.close()

######################################################################################
######################## SAMPLE 2 BAR PLOT ###########################################
######################################################################################

# create bar plot
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, figsize=(20,6), dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(df2_amplicon))
width = 0.2

bar1 = ax.bar(x_axis, df2_amplicon['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2 = ax.bar(x_axis + width, df2_amplicon['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, df2_amplicon['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, df2_amplicon['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, df2_amplicon['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, df2_amplicon['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax2.set_ylim(0, 10)
ax.set_ylim(10, 100)
ax2.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)

d = .005  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_ylabel('Counts', loc='top')
ax2.set_xticks(x_axis + width,df2_amplicon['start'], fontsize=8, rotation=90)
ax.set_title(sample2.split("_NT",1)[0])
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )
plt.savefig(sample2 + '_barplot_non_canonical.jpeg')
plt.close()
# plt.show()

# df2_amplicon.to_csv(sample2 + "_novel_sgRNA_softwares_table.csv", sep=';', header=True, index=False)

######################################################################################
######################## SAMPLE 2 CONCORDANCE RATE ORFS REGIONS ######################
######################################################################################
df2_ORFs_starting = pd.read_csv(r'intersection_ORFS_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end", "ORF"])

merge2 = df2_ORFs_starting.merge(df2_amplicon, how='outer')
conte2 = merge2[['ORF','counts_periscope','counts_leTRS','sgDI_count']].groupby('ORF').sum().reset_index()

df2_ORFs = df_ORFs_names.merge(conte2, on='ORF', how='outer').fillna(0)
   
df2_amplicon_ORFS = df2_amplicon.merge(df_ORFs_starting, how='outer')
df2_amplicon_ORFS.to_csv(sample2 + "_novel_sgRNA_softwares_table.csv", sep=';', header=True, index=False)

######################################################################################
######################## SAMPLE 2 BAR PLOT ORFS REGIONS ##############################
######################################################################################

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(df2_ORFs))
width = 0.2

bar1 = ax.bar(x_axis, df2_ORFs['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2= ax.bar(x_axis + width, df2_ORFs['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, df2_ORFs['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, df2_ORFs['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, df2_ORFs['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, df2_ORFs['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax2.set_ylim(0, 10)
ax.set_ylim(10, 100)
ax2.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xticks(x_axis+width,df2_ORFs['ORF'], fontsize=7, rotation=90)
ax2.set_ylabel('Counts', loc='top')
ax.set_title(sample2.split("_NT",1)[0])
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )
# plt.show()
plt.savefig(sample2 + '_barplot_non_canonical_ORFS.jpeg')
plt.close()







######################################################################################
############################ MEAN BETWEEN REPLICATES #################################
######################################################################################
sub_df = df_amplicon[['chromosome','start','end','counts_periscope','counts_leTRS','sgDI_count']]
sub_df2 = df2_amplicon[['chromosome','start','end','counts_periscope','counts_leTRS','sgDI_count']]
mergione_df = pd.concat((sub_df, sub_df2))
mergione_df  = mergione_df .groupby(['chromosome','start','end'], sort=False).mean()
mergione_df  = mergione_df .reset_index()

name = sample.split("_1",1)[0] + '_' + sample.split("_",3)[2]

# create bar plot
fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, figsize=(20,6), dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(mergione_df))
width = 0.2

bar1 = ax.bar(x_axis, mergione_df['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2 = ax.bar(x_axis + width, mergione_df['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, mergione_df['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, mergione_df['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, mergione_df['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, mergione_df['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax2.set_ylim(0, 10)
ax.set_ylim(10, 100)
ax2.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)

d = .005  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_ylabel('Counts', loc='top')
ax2.set_xticks(x_axis + width,mergione_df['start'], fontsize=8, rotation=90)
ax.set_title(name)
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )
plt.savefig(name + '_barplot_non_canonical_means.jpeg')
plt.close()
# plt.show()



########################## ORFS ################################

sub_df_ORFs = df_ORFs[['ORF','counts_periscope','counts_leTRS','sgDI_count']]
sub_df2_ORFs = df2_ORFs[['ORF','counts_periscope','counts_leTRS','sgDI_count']]
merge_sub_ORFs = pd.concat((sub_df_ORFs, sub_df2_ORFs))
merge_sub_ORFs = merge_sub_ORFs.groupby('ORF', sort=False).mean()
merge_sub_ORFs = merge_sub_ORFs.reset_index()


fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(merge_sub_ORFs))
width = 0.2

bar1 = ax.bar(x_axis, merge_sub_ORFs['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar2= ax.bar(x_axis + width, merge_sub_ORFs['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar3 = ax.bar(x_axis + width*2, merge_sub_ORFs['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

bar4 = ax2.bar(x_axis, merge_sub_ORFs['counts_periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
bar5 = ax2.bar(x_axis + width, merge_sub_ORFs['counts_leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
bar6 = ax2.bar(x_axis + width*2, merge_sub_ORFs['sgDI_count'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

ax2.set_ylim(0, 10)
ax.set_ylim(10, 100)
ax2.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False, top=False)  # don't put tick labels at the top
ax2.tick_params(labeltop=False, bottom=False)

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

ax2.set_xticks(x_axis+width,merge_sub_ORFs['ORF'], fontsize=7, rotation=90)
ax2.set_ylabel('Counts', loc='top')
ax.set_title(name)
ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )
# plt.show()
plt.savefig(name + '_barplot_non_canonical_ORFS.jpeg')
plt.close()

# save means to use for orf plot
merge_sub_ORFs.to_csv(name +'_mean_orfs_nc-sgRNA.csv', index = False, header=True, sep = "\t")
