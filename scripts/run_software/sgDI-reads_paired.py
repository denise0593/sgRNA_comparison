#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:46:50 2022

@author: 'Denise Lavezzari'
"""

## The script takes as input the results from the sgDI-tector from the two reads and merge the results.

import os
import pandas as pd
import sys


os.chdir(os.getcwd())
sample = sys.argv[1]

os.chdir(sample)

sgDI_R1 = pd.read_csv("sgDI_results/R1/sgRNA_junctions.csv")
sgDI_R2 = pd.read_csv("sgDI_results/R2/sgRNA_junctions.csv")

merged = pd.concat([sgDI_R1,sgDI_R2]).sort_values('counts').drop_duplicates('ri_pos', keep='last')

orf = merged["ORF"].str.contains("ORF")
orf_df = merged[orf]

nc = merged[~orf]

orf_df.to_csv("sgDI_results/canonical_sgDI.csv", sep=',', header=True, index=False)
nc.to_csv("sgDI_results/non_canonical_sgDI.csv", sep=',', header=True, index=False)
