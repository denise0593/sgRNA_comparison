#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:00:46 2022

@author: 'Denise Lavezzari'
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:17:52 2022

@author: 'Denise Lavezzari'
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted
import seaborn as sns

os.chdir("NT")

sample = "6264_1_T0_NT_non_trattate"
sample2 = "6264_2_T0_NT_non_trattate"

# os.chdir(os.getcwd())
# sample = sys.argv[1]
# sample2 = sys.argv[2]

# orf dataframe
orfs = pd.DataFrame()
df_bed = pd.DataFrame()

name = sample.split("_1",1)[0] + '_' + sample.split("_",3)[2]

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
############################### SAMPLE 1 COUNTS ######################################
######################################################################################
df_amplicon_per = pd.read_csv(r'softwares_merged_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end"])
df_amplicon_leTRS = df_amplicon_per.copy()
df_amplicon_sgDI = df_amplicon_per.copy()


df_amplicon_per[['Counts']] = ""
df_amplicon_leTRS[['Counts']] = ""
df_amplicon_sgDI[['Counts']] = ""

set_periscope = set()
set_leTRS = set()
set_sgDI = set()

# count sgRNA contained in the intervals, counts of reads supporting the intervals
for i in range(len(df_amplicon_per)):
    start = df_amplicon_per.loc[i,'start']
    end = df_amplicon_per.loc[i,'end']
    if sample_periscope['pos'].between(start,end).any():
        df_amplicon_per.loc[i,'Counts'] = sample_periscope[sample_periscope['pos'].between(start,end)]['sgRNA_count'].sum()
        set_periscope.add(df_amplicon_per.loc[i,'start'])
    else: 
        df_amplicon_per.loc[i,'Counts'] = 0
        
    if sample_leTRS['pos'].between(start,end).any():
        df_amplicon_leTRS.loc[i,'Counts'] = sample_leTRS[sample_leTRS['pos'].between(start,end)]['nb_count'].sum()
        set_leTRS.add(df_amplicon_leTRS.loc[i,'start'])
    else: 
        df_amplicon_leTRS.loc[i,'Counts'] = 0
        
    if sample_sgDI['pos'].between(start,end).any():
        df_amplicon_sgDI.loc[i,'Counts'] = sample_sgDI[sample_sgDI['pos'].between(start,end)]['counts'].sum()
        set_sgDI.add(df_amplicon_sgDI.loc[i,'start'])
    else: 
        df_amplicon_sgDI.loc[i,'Counts'] = 0
        
df_amplicon_per["Software"] = "Periscope"
df_amplicon_leTRS["Software"] = "LeTRS" 
df_amplicon_sgDI["Software"] = "sgDI-tector"


######################################################################################
############################### SAMPLE 2 COUNTS ######################################
######################################################################################
df2_amplicon_per = pd.read_csv(r'softwares_merged_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end"])

df2_amplicon_leTRS = df2_amplicon_per.copy()
df2_amplicon_sgDI = df2_amplicon_per.copy()

df2_amplicon_per[['Counts']] = ""
df2_amplicon_leTRS[['Counts']] = ""
df2_amplicon_sgDI[['Counts']] = ""

set2_periscope = set()
set2_leTRS = set()
set2_sgDI = set()

# count sgRNA contained in the intervals, counts of reads supporting the intervals
for i in range(len(df2_amplicon_per)):
    start = df2_amplicon_per.loc[i,'start']
    end = df2_amplicon_per.loc[i,'end']
    if sample2_periscope['pos'].between(start,end).any():
        df2_amplicon_per.loc[i,'Counts'] = sample2_periscope[sample2_periscope['pos'].between(start,end)]['sgRNA_count'].sum()
        set2_periscope.add(df2_amplicon_per.loc[i,'start'])
    else: 
        df2_amplicon_per.loc[i,'Counts'] = 0
        
    if sample2_leTRS['pos'].between(start,end).any():
        df2_amplicon_leTRS.loc[i,'Counts'] = sample2_leTRS[sample2_leTRS['pos'].between(start,end)]['nb_count'].sum()
        set2_leTRS.add(df2_amplicon_leTRS.loc[i,'start'])
    else: 
        df2_amplicon_leTRS.loc[i,'Counts'] = 0
        
    if sample2_sgDI['pos'].between(start,end).any():
        df2_amplicon_sgDI.loc[i,'Counts'] = sample2_sgDI[sample2_sgDI['pos'].between(start,end)]['counts'].sum()
        set2_sgDI.add(df2_amplicon_sgDI.loc[i,'start'])
    else: 
        df2_amplicon_sgDI.loc[i,'Counts'] = 0
        
df2_amplicon_per["Software"] = "Periscope"
df2_amplicon_leTRS["Software"] = "LeTRS" 
df2_amplicon_sgDI["Software"] = "sgDI-tector"

full_df=pd.concat((df_amplicon_per,df_amplicon_leTRS,df_amplicon_sgDI,df2_amplicon_per,df2_amplicon_leTRS,df2_amplicon_sgDI))



############################### PLOT ######################################

fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4,6]}, figsize=(20,8), dpi=600)

fig.subplots_adjust(hspace=0.05)

x_axis = np.arange(len(full_df['start'].unique()))
width = 0.2

sns.barplot(x="start", y="Counts",  hue="Software", data=full_df, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black', ax=axes[0])
sns.swarmplot(x="start", y="Counts", hue="Software", data=full_df, palette='dark:0', dodge=True, edgecolor="gray", legend=False, ax=axes[0])

sns.barplot(x="start", y="Counts",  hue="Software", data=full_df, palette=["darkgreen","blue","red"], alpha = 0.4,  edgecolor = 'black', ax=axes[1])
sns.swarmplot(x="start", y="Counts", hue="Software", data=full_df, palette='dark:0', dodge=True,  edgecolor="gray", legend=False, ax=axes[1])


# remove legend repetition
axes[1].get_legend().set_visible(False)

# remove y label repetition an set y label
axes[1].set(ylabel=None) 
axes[0].set(ylabel=None) 
axes[1].set_ylabel('Counts', loc='top', fontsize=15)

# remove x label repetition and set the one on the bottom subplot
axes[0].set(xlabel=None) 
axes[1].set(xlabel=None) 
axes[0].set_xticks([])
axes[1].set_xticklabels(full_df['start'].unique(), rotation = 90, fontsize=7)

axes[1].set_ylim(0, 10)
axes[0].set_ylim(10,100)  
 
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

fig.suptitle(name, fontsize=15)

plt.savefig(name + '_plot_counts_non_canonical.jpeg')
plt.close()



######################################################################################
######################## SAMPLE 1 CONCORDANCE RATE ORFS REGIONS ######################
######################################################################################
df_ORFs_starting = pd.read_csv(r'intersection_ORFS_'+ sample + '_' + sample2 + '.bed', delimiter='\t', names=["chromosome", "start","end", "ORF"])

rep1_df=pd.concat((df_amplicon_per,df_amplicon_leTRS,df_amplicon_sgDI))
rep2_df=pd.concat((df2_amplicon_per,df2_amplicon_leTRS,df2_amplicon_sgDI))

merge_rep1 = df_ORFs_starting.merge(rep1_df, how='outer')
merge_rep2 = df_ORFs_starting.merge(rep2_df, how='outer')

merge_rep1 = merge_rep1[['ORF', 'Counts', 'Software']].groupby(['ORF','Software']).sum().reset_index()
merge_rep2 = merge_rep2[['ORF', 'Counts', 'Software']].groupby(['ORF','Software']).sum().reset_index()


df_ORFs_names = pd.DataFrame()
df_ORFs_names['ORF'] = ['ORF 1ab','ORF 1ab','ORF 1ab','ORF S','ORF S','ORF S','ORF 3a','ORF 3a','ORF 3a','ORF E','ORF E','ORF E','ORF M','ORF M','ORF M','ORF 6','ORF 6', 'ORF 6','ORF 7a','ORF 7a','ORF 7a', 'ORF 7b', 'ORF 7b', 'ORF 7b', 'ORF 8', 'ORF 8', 'ORF 8', 'ORF10','ORF10','ORF10','ORF N','ORF N', 'ORF N']
df_ORFs_names['Software'] = ['Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector','Periscope','LeTRS','sgDI-tector']


df_ORFs_rep1 = df_ORFs_names.merge(merge_rep1, on=['ORF','Software'], how='outer').fillna(0)
df_ORFs_rep2 = df_ORFs_names.merge(merge_rep2, on=['ORF','Software'], how='outer').fillna(0)

# add to dataframe row with no counts for any software
# for i in df_ORFs_rep1.index:
#     if df_ORFs_rep1['Software'][i] == 0:
#         new_row_per = pd.DataFrame({'ORF': df_ORFs_rep1['ORF'][i], 'Software':'Periscope', 'Counts':[0]})
#         new_row_LeTRS = pd.DataFrame({'ORF':df_ORFs_rep1['ORF'][i], 'Software':'LeTRS', 'Counts':[0]})
#         new_row_sgDI = pd.DataFrame( {'ORF':df_ORFs_rep1['ORF'][i], 'Software':'sgDI-tector', 'Counts':[0]})        
#         df_ORFs_rep1 = pd.concat((df_ORFs_rep1, new_row_per))
#         df_ORFs_rep1 = pd.concat((df_ORFs_rep1, new_row_LeTRS))
#         df_ORFs_rep1 = pd.concat((df_ORFs_rep1, new_row_sgDI))
        
# df_ORFs_rep1 = df_ORFs_rep1[df_ORFs_rep1.Software != 0]

# # add to dataframe row with no counts for any software
# for i in df_ORFs_rep2.index:
#     if df_ORFs_rep2['Software'][i] == 0:
#         new_row_per = pd.DataFrame({'ORF': df_ORFs_rep2['ORF'][i], 'Software':'Periscope', 'Counts':[0]})
#         new_row_LeTRS = pd.DataFrame({'ORF':df_ORFs_rep2['ORF'][i], 'Software':'LeTRS', 'Counts':[0]})
#         new_row_sgDI = pd.DataFrame( {'ORF':df_ORFs_rep2['ORF'][i], 'Software':'sgDI-tector', 'Counts':[0]})        
#         df_ORFs_rep2 = pd.concat((df_ORFs_rep2, new_row_per))
#         df_ORFs_rep2 = pd.concat((df_ORFs_rep2, new_row_LeTRS))
#         df_ORFs_rep2 = pd.concat((df_ORFs_rep2, new_row_sgDI))
        
# df_ORFs_rep2 = df_ORFs_rep2[df_ORFs_rep2.Software != 0]

df_ORFs_full = pd.concat((df_ORFs_rep1,df_ORFs_rep2))

order = ['Periscope','LeTRS','sgDI-tector']

df = df_ORFs_full.sort_values([])

############################### PLOT ######################################

fig.subplots_adjust(hspace=0.05) 
fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4,6]}, dpi=600)
fig.subplots_adjust(hspace=0.05) 

x_axis = np.arange(len(df_ORFs_full['ORF'].unique()))
width = 0.2

sns.barplot(x="ORF", y="Counts",  hue="Software", data=df_ORFs_full, palette=["darkgreen","blue","red"], alpha = 0.4, edgecolor = 'black', ax=axes[0])
sns.swarmplot(x="ORF", y="Counts", hue="Software", data=df_ORFs_full, palette='dark:0', dodge=True, edgecolor="gray", legend=False, ax=axes[0], size = 4)

sns.barplot(x="ORF", y="Counts",  hue="Software", data=df_ORFs_full, palette=["darkgreen","blue","red"], alpha = 0.4, edgecolor = 'black', ax=axes[1])
sns.swarmplot(x="ORF", y="Counts", hue="Software", data=df_ORFs_full, palette='dark:0', dodge=True,  edgecolor="gray", legend=False, ax=axes[1], size = 4)


# remove legend repetition
axes[1].get_legend().set_visible(False)

# remove y label repetition an set y label
axes[1].set(ylabel=None) 
axes[0].set(ylabel=None) 
axes[1].set_ylabel('Counts', loc='top')

# remove x label repetition and set the one on the bottom subplot
axes[0].set(xlabel=None) 
axes[1].set(xlabel=None) 
axes[0].set_xticks([])
axes[1].set_xticklabels(df_ORFs_full['ORF'].unique(), rotation = 90, fontsize=7)

axes[1].set_ylim(0, 10)
axes[0].set_ylim(10,100)  
 
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

fig.suptitle(name)

plt.savefig(name + '_plot_counts_non_canonical_ORFs.jpeg')
plt.close()
