#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:22:53 2023

@author: 'Denise Lavezzari'
"""

import pandas as pd
import os
import plot_orf_in_time as pot

# THE SCRIPT TAKES AS INPUT THE CSV FILES OF EACH TIME POINT AND CREATES AS OUTPUT A PLOT FOR EACH ORFS AT DIFFERENT TIME POINTS

# folder with csv files
os.chdir("Files")

def listaDati(pos,string):
    lista = [df_0h.loc[pos,string],df_24h.loc[pos,string],df_48h.loc[pos,string],df_72h.loc[pos,string],df_96h.loc[pos,string]]
    return lista

###################################
############ c-sgRNAs #############
###################################

c_df_0h = pd.read_csv('6264_mean_orfs.csv', sep=",")
c_df_24h = pd.read_csv('6265_mean_orfs.csv', sep=",")
c_df_48h = pd.read_csv('6266_mean_orfs.csv', sep=",")
c_df_72h = pd.read_csv('6267_mean_orfs.csv', sep=",")
c_df_96h = pd.read_csv('6268_mean_orfs.csv', sep=",")

sample_column = ['6264_0h','6265_24h','6266_48h','6267_72h','6268_96h']


# ORF 1ab
c_df_ORF1ab = pd.DataFrame()
c_df_ORF1ab['ORF1 ab'] = sample_column
c_df_ORF1ab['Periscope'] = listaDati(0,'counts_periscope')
c_df_ORF1ab['leTRS'] = listaDati(0,'counts_leTRS')
c_df_ORF1ab['sgDI'] = listaDati(0,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF1ab)

# ORF S
c_df_ORFS = pd.DataFrame()
c_df_ORFS['ORF S'] = sample_column
c_df_ORFS['Periscope'] = listaDati(1,'counts_periscope')
c_df_ORFS['leTRS'] = listaDati(1,'counts_leTRS')
c_df_ORFS['sgDI'] = listaDati(1,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORFS)

# ORF 3a
c_df_ORF3a = pd.DataFrame()
c_df_ORF3a['ORF 3a'] = sample_column
c_df_ORF3a['Periscope'] = listaDati(2,'counts_periscope')
c_df_ORF3a['leTRS'] = listaDati(2,'counts_leTRS')
c_df_ORF3a['sgDI'] = listaDati(2,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF3a)

# ORF E
c_df_ORFE = pd.DataFrame()
c_df_ORFE['ORF E'] = sample_column
c_df_ORFE['Periscope'] = listaDati(3,'counts_periscope')
c_df_ORFE['leTRS'] = listaDati(3,'counts_leTRS')
c_df_ORFE['sgDI'] = listaDati(3,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORFE)

# ORF M
c_df_ORFM = pd.DataFrame()
c_df_ORFM['ORF M'] = sample_column
c_df_ORFM['Periscope'] = listaDati(4,'counts_periscope')
c_df_ORFM['leTRS'] = listaDati(4,'counts_leTRS')
c_df_ORFM['sgDI'] = listaDati(4,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORFM)

# ORF 6
c_df_ORF6 = pd.DataFrame()
c_df_ORF6['ORF 6'] = sample_column
c_df_ORF6['Periscope'] = listaDati(5,'counts_periscope')
c_df_ORF6['leTRS'] = listaDati(5,'counts_leTRS')
c_df_ORF6['sgDI'] = listaDati(5,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF6)

# ORF 7a
c_df_ORF7a = pd.DataFrame()
c_df_ORF7a['ORF 7a'] = sample_column
c_df_ORF7a['Periscope'] = listaDati(6,'counts_periscope')
c_df_ORF7a['leTRS'] = listaDati(6,'counts_leTRS')
c_df_ORF7a['sgDI'] = listaDati(6,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF7a)

# ORF 7b
c_df_ORF7b = pd.DataFrame()
c_df_ORF7b['ORF 7b'] = sample_column
c_df_ORF7b['Periscope'] = listaDati(7,'counts_periscope')
c_df_ORF7b['leTRS'] = listaDati(7,'counts_leTRS')
c_df_ORF7b['sgDI'] = listaDati(7,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF7b)

# ORF 8
c_df_ORF8 = pd.DataFrame()
c_df_ORF8['ORF 8'] = sample_column
c_df_ORF8['Periscope'] = listaDati(8,'counts_periscope')
c_df_ORF8['leTRS'] = listaDati(8,'counts_leTRS')
c_df_ORF8['sgDI'] = listaDati(8,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF8)

# ORF 10
c_df_ORF10 = pd.DataFrame()
c_df_ORF10['ORF 10'] = sample_column
c_df_ORF10['Periscope'] = listaDati(9,'counts_periscope')
c_df_ORF10['leTRS'] = listaDati(9,'counts_leTRS')
c_df_ORF10['sgDI'] = listaDati(9,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORF10)

# ORF N
c_df_ORFN = pd.DataFrame()
c_df_ORFN['ORF N'] = sample_column
c_df_ORFN['Periscope'] = listaDati(10,'counts_periscope')
c_df_ORFN['leTRS'] = listaDati(10,'counts_leTRS')
c_df_ORFN['sgDI'] = listaDati(10,'sgDI_count')

pot.create_ORF_time_csgRNA(c_df_ORFN)


###################################
############ nc-sgRNAs ############
###################################

# open csv file with mean values for each orfs for each time point
df_0h = pd.read_csv('6264_T0_mean_orfs_nc-sgRNA.csv', sep="\t")
df_24h = pd.read_csv('6265_24h_mean_orfs_nc-sgRNA.csv', sep="\t")
df_48h = pd.read_csv('6266_48h_mean_orfs_nc-sgRNA.csv', sep="\t")
df_72h = pd.read_csv('6267_72h_mean_orfs_nc-sgRNA.csv', sep="\t")
df_96h = pd.read_csv('6268_96h_mean_orfs_nc-sgRNA.csv', sep="\t")

sample_column = ['6264_0h','6265_24h','6266_48h','6267_72h','6268_96h']

# ORF 1ab
df_ORF1ab = pd.DataFrame()
df_ORF1ab['ORF1 ab'] = sample_column
df_ORF1ab['Periscope'] = listaDati(0,'counts_periscope')
df_ORF1ab['leTRS'] = listaDati(0,'counts_leTRS')
df_ORF1ab['sgDI'] = listaDati(0,'sgDI_count')

# Recall the script plot_orf_in_time to create the plot
pot.create_ORF_time_ncsgRNA(df_ORF1ab)

# ORF S
df_ORFS = pd.DataFrame()
df_ORFS['ORF S'] = sample_column
df_ORFS['Periscope'] = listaDati(1,'counts_periscope')
df_ORFS['leTRS'] = listaDati(1,'counts_leTRS')
df_ORFS['sgDI'] = listaDati(1,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORFS)

# ORF 3a
df_ORF3a = pd.DataFrame()
df_ORF3a['ORF 3a'] = sample_column
df_ORF3a['Periscope'] = listaDati(2,'counts_periscope')
df_ORF3a['leTRS'] = listaDati(2,'counts_leTRS')
df_ORF3a['sgDI'] = listaDati(2,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORF3a)

# ORF E
df_ORFE = pd.DataFrame()
df_ORFE['ORF E'] = sample_column
df_ORFE['Periscope'] = listaDati(3,'counts_periscope')
df_ORFE['leTRS'] = listaDati(3,'counts_leTRS')
df_ORFE['sgDI'] = listaDati(3,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORFE)

# ORF M
df_ORFM = pd.DataFrame()
df_ORFM['ORF M'] = sample_column
df_ORFM['Periscope'] = listaDati(4,'counts_periscope')
df_ORFM['leTRS'] = listaDati(4,'counts_leTRS')
df_ORFM['sgDI'] = listaDati(4,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORFM)

# ORF 6
df_ORF6 = pd.DataFrame()
df_ORF6['ORF 6'] = sample_column
df_ORF6['Periscope'] = listaDati(5,'counts_periscope')
df_ORF6['leTRS'] = listaDati(5,'counts_leTRS')
df_ORF6['sgDI'] = listaDati(5,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORF6)

# ORF 7a
df_ORF7a = pd.DataFrame()
df_ORF7a['ORF 7a'] = sample_column
df_ORF7a['Periscope'] = listaDati(6,'counts_periscope')
df_ORF7a['leTRS'] = listaDati(6,'counts_leTRS')
df_ORF7a['sgDI'] = listaDati(6,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORF7a)

# ORF 7b
df_ORF7b = pd.DataFrame()
df_ORF7b['ORF 7b'] = sample_column
df_ORF7b['Periscope'] = listaDati(7,'counts_periscope')
df_ORF7b['leTRS'] = listaDati(7,'counts_leTRS')
df_ORF7b['sgDI'] = listaDati(7,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORF7b)

# ORF 8
df_ORF8 = pd.DataFrame()
df_ORF8['ORF 8'] = sample_column
df_ORF8['Periscope'] = listaDati(8,'counts_periscope')
df_ORF8['leTRS'] = listaDati(8,'counts_leTRS')
df_ORF8['sgDI'] = listaDati(8,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORF8)

# ORF 10
df_ORF10 = pd.DataFrame()
df_ORF10['ORF 10'] = sample_column
df_ORF10['Periscope'] = listaDati(9,'counts_periscope')
df_ORF10['leTRS'] = listaDati(9,'counts_leTRS')
df_ORF10['sgDI'] = listaDati(9,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORF10)

# ORF N
df_ORFN = pd.DataFrame()
df_ORFN['ORF N'] = sample_column
df_ORFN['Periscope'] = listaDati(10,'counts_periscope')
df_ORFN['leTRS'] = listaDati(10,'counts_leTRS')
df_ORFN['sgDI'] = listaDati(10,'sgDI_count')

pot.create_ORF_time_ncsgRNA(df_ORFN)

