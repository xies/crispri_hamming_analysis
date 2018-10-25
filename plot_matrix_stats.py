#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 18:44:12 2018

@author: xies
"""

import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from Bio import SeqIO
from itertools import count, izip

filename = '/data/crispri_hamming/big_distance_matrix.csv'
Dists = np.zeros(19172)
Dists = np.loadtxt(filename,delimiter=',')

D = Dists[:,0]
index = 0

#with open(filename, 'rb') as f:
#    lines = f.readlines()
#    for l in lines:
#        myarray = np.fromstring(l, dtype=float, sep=',')
#        
#        D[index] = myarray[index + 1:].min()
#        
#        index += 1
        
#min_distances = D.min()

plt.figure()
bins = range( int(D.max()) +2)
plt.hist(D,bins=bins, normed=True, cumulative=False, histtype='bar',align='left')
plt.xlabel('Minimum mismatches')
plt.ylabel('CDF')


# Check manually for exact matches
filename = '/data/ForMimi/AllSgRNAsOct4'
seqs = [rec for rec in SeqIO.parse(filename,'fasta')]

zeros = np.where(D == 0)[0]
matches = Dists[D == 0,1].astype(np.int)

df = []
for (i, matching_pair) in izip(count(), izip(zeros,matches)):
    df.append((i,matching_pair[0],matching_pair[1],
               str(seqs[matching_pair[0]].seq), str(seqs[matching_pair[1]].seq) ))
df = pd.DataFrame(df)
df.to_csv('/data/crispri_hamming/exact_matches.txt',sep='\t')
