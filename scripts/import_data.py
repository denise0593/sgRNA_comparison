#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:37:11 2023

@author: 'Denise Lavezzari'
"""

import pandas as pd

# The script define the function to import the data from the output of the different software and standardize the orf name

def periscope_import(sample):
    # periscope
    sample_periscope = pd.read_csv(sample + '_periscope_counts.csv')
    sample_periscope['orf'] = sample_periscope['orf'].str.replace(r'ORF', '')
    sample_periscope['orf'] = 'ORF ' + sample_periscope['orf'].str.upper()
    sample_periscope['orf'] = sample_periscope['orf'].replace({'ORF 1A': 'ORF 1ab'})
    sample_periscope['orf'] = sample_periscope['orf'].replace({'ORF 3A': 'ORF 3a'})
    sample_periscope['orf'] = sample_periscope['orf'].replace({'ORF 7A': 'ORF 7a'})
    sample_periscope['orf'] = sample_periscope['orf'].replace({'ORF 7B': 'ORF 7b'})
    return sample_periscope

def leTRS_import(sample):
    # leTRS
    sample_leTRS = pd.read_csv('LeTRS_output_cov0/results/known_junction.tab', sep="\t")
    sample_leTRS.drop(sample_leTRS.tail(7).index, inplace=True)
    sample_leTRS['subgenome'] = sample_leTRS['subgenome'].str.replace(r'ORF', '')
    sample_leTRS['subgenome'] = 'ORF ' + sample_leTRS['subgenome'].str.upper()
    sample_leTRS['subgenome'] = sample_leTRS['subgenome'].replace({'ORF 1AB': 'ORF 1ab'})
    sample_leTRS['subgenome'] = sample_leTRS['subgenome'].replace({'ORF 3A': 'ORF 3a'})
    sample_leTRS['subgenome'] = sample_leTRS['subgenome'].replace({'ORF 7A': 'ORF 7a'})
    sample_leTRS['subgenome'] = sample_leTRS['subgenome'].replace({'ORF 7B': 'ORF 7b'})
    sample_leTRS['peak_count'] = sample_leTRS['peak_count'].str.split("(", 1).str.get(0).astype(str).astype(int)
    sample_leTRS = sample_leTRS[ sample_leTRS.peak_count != 0 ]
    return sample_leTRS

def sgDI_import(sample):
    # sgDI-tector
    sample_sgDI = pd.read_csv('sgDI_results/canonical_sgDI.csv', sep=",")
    sample_sgDI['ORF'] = sample_sgDI['ORF'].str.upper()
    sample_sgDI['ORF'] = sample_sgDI['ORF'].replace({'ORF 1AB': 'ORF 1ab'})
    sample_sgDI['ORF'] = sample_sgDI['ORF'].replace({'ORF 3A': 'ORF 3a'})
    sample_sgDI['ORF'] = sample_sgDI['ORF'].replace({'ORF 7A': 'ORF 7a'})
    sample_sgDI['ORF'] = sample_sgDI['ORF'].replace({'ORF 7B': 'ORF 7b'})
    return sample_sgDI

def mergione(df_periscope,df_leTRS,df_sgDI):
    # create a subset of data
    sub_periscope = pd.DataFrame()
    sub_periscope[['orf', 'counts']] = df_periscope[['orf', 'sgRNA_count']]
    sub_leTRS = pd.DataFrame()
    sub_leTRS[['orf', 'counts']] = df_leTRS[['subgenome', 'peak_count']]
    sub_sgDI = pd.DataFrame()
    sub_sgDI[['orf', 'counts']] = df_sgDI[['ORF', 'counts']]
    orfs = pd.DataFrame()
    orfs['orf'] = ['ORF 1ab','ORF S','ORF 3a','ORF E','ORF M','ORF 6', 'ORF 7a', 'ORF 7b', 'ORF 8', 'ORF 10']
    mergione = orfs.merge(sub_periscope, on='orf', how='outer').merge(sub_leTRS,on='orf', how='outer', suffixes=('_periscope', '_leTRS')).merge(sub_sgDI,on='orf', how='outer').fillna(0)
    return mergione
    
