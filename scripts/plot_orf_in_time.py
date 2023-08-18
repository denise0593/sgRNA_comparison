#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:22:21 2023

@author: 'Denise Lavezzari'
"""

# import os
# import sys
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def create_ORF_time_csgRNA(orf_df): 
    name = orf_df.columns[0]
    
    # create bar plot
    fig, (ax, ax2, ax3) = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [4, 3, 3]}, dpi=600)
    fig.subplots_adjust(hspace=0.05) 
    
    x_axis = np.arange(len(orf_df))
    width = 0.2
    
    # define axes for plot
    bar1 = ax.bar(x_axis, orf_df['Periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
    bar2 = ax.bar(x_axis + width, orf_df['leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
    bar3 = ax.bar(x_axis + width*2, orf_df['sgDI'], width, edgecolor = 'black', color = 'red', alpha = 0.4)
        
    bar4 = ax2.bar(x_axis, orf_df['Periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
    bar5 = ax2.bar(x_axis + width, orf_df['leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
    bar6 = ax2.bar(x_axis + width*2, orf_df['sgDI'], width, edgecolor = 'black', color = 'red', alpha = 0.4)
    
    bar7 = ax3.bar(x_axis, orf_df['Periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
    bar8 = ax3.bar(x_axis + width, orf_df['leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
    bar9 = ax3.bar(x_axis + width*2, orf_df['sgDI'], width, edgecolor = 'black', color = 'red', alpha = 0.4)
    
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
    
    ax.set_xticks(x_axis+width,orf_df[name], fontsize=7)
    ax.set_ylabel('Counts', loc='center')
    ax.set_title(name)
    ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector'), loc='upper left')
        
    # plt.show()
    plt.savefig(name + '_barplot_in_time.jpeg')
    plt.close()
    

def create_ORF_time_ncsgRNA(orf_df): 
    name = orf_df.columns[0]
    
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,6]}, dpi=600)
    fig.subplots_adjust(hspace=0.05) 

    x_axis = np.arange(len(orf_df))
    width = 0.2

    bar1 = ax.bar(x_axis, orf_df['Periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
    bar2= ax.bar(x_axis + width, orf_df['leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
    bar3 = ax.bar(x_axis + width*2, orf_df['sgDI'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

    bar4 = ax2.bar(x_axis, orf_df['Periscope'], width, edgecolor = 'black', color = 'darkgreen', alpha = 0.4)
    bar5 = ax2.bar(x_axis + width, orf_df['leTRS'], width, edgecolor = 'black', color = 'blue', alpha = 0.4)
    bar6 = ax2.bar(x_axis + width*2, orf_df['sgDI'], width, edgecolor = 'black', color = 'red', alpha = 0.4)

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

    ax2.set_xticks(x_axis+width,orf_df[name], fontsize=7)
    ax2.set_ylabel('Counts', loc='top')
    ax.set_title(name)
    ax.legend( (bar1, bar2, bar3), ('Periscope', 'LeTRS', 'sgDI-tector') )

    # plt.show()
    plt.savefig(str(name) + '_barplot_in_time_nc-sgRNA.jpeg')
    plt.close()
