#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 20:43:45 2018

@author: xies
"""

import numpy as np
import matplotlib.pylab as plt
from Bio import SeqIO, SeqRecord
from itertools import count, izip

guide_window = 75

# Load ChIP peaks
filename = '/data/crispri_hamming/oct4_chip_flanking.fa'
peaks = [rec for rec in SeqIO.parse(filename,'fasta')]
Npeaks = len(peaks)

# Merge ChIP peaks that are within 75bp


# Load sgRNA seqs
filename = '/data/crispri_hamming/sgRNA_unique.fasta'
guides = [rec for rec in SeqIO.parse(filename,'fasta')]
Nguides = len(guides)
# Fix the sequence name
for g in guides:
    g.name = get_peak_name(g.description)

# Load closest hamming distance
filename = '/data/crispri_hamming/big_distance_matrix.csv'
D = np.loadtxt(filename,delimiter=',')

# Go through all peaks and find the 'most unique' guide
guide_names = [g.name for g in guides]
max_distances = np.ones(Npeaks)*-1
for (i,p) in enumerate(peaks):
    # Find all sgRNA belonging to this peak
    indices = [j for j,n in enumerate(guide_names) if n == p.name]
    if len(indices) > 0:
        max_distances[i] = D[indices,0].max()
    
# Plot the inverse distribution of whether a ChIP peak would have a sgRNA "unique" enough
outs = plt.hist(max_distances[max_distances >= 0],cumulative=True,normed=True)

def peaks_overlapping(p1,p2,w):
    start1 = 
    return

def get_peak_chromosome(name):
    return name.split(':')[0]

def get_peak_location(name):
    start = name.split(':')[1]
    end = start.split('-')[1]
    start = start.split('-')[0]
    return start, end

def get_peak_name(desc):
    return desc.split(' ')[1]

