#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 17:30:42 2018

@author: xies
"""

import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from Bio import SeqIO, motifs, SeqRecord
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from timeit import default_timer

filename = '/data/crispri_hamming/nanog_chip_peaks.fa'
peaks = [rec for rec in SeqIO.parse(filename,'fasta')]
for rec in peaks:
    rec.seq.alphabet = IUPACUnambiguousDNA()
Npeaks = len(peaks)
has_guide = np.zeros(Npeaks,dtype=bool)
Lpeak = len(peaks[0])

# Load NANOG PWM and rewrite into JASPAR format
filename = '/data/crispri_hamming/nanog_GSE11724.jaspar'
with open(filename)as fh:
    m = [m for m in motifs.parse(fh, "jaspar")]
pwm = m[0]
pssm = pwm.pssm

Lmotif = len(pwm)

# Assuming 25-25-25-25 background, find min and max scores and 80% threshold
min_score = pssm.min; max_score = pssm.max
threshold_80perc = (max_score - min_score) * 0.80 + min_score

# Search through sequences and find hits to PSSM
# Then pull sequences with flanking arms
flanking_size = 30

tic = default_timer()
sites = []
for ts in peaks:
    # DO NOT use 'both=True' because negative indexing on pssm.search is broken AF
    chr_name = get_peak_chromosome(ts.name)
    peak_start = get_peak_location(ts.name)[0]
    peak_end = get_peak_location(ts.name)[1]
    
    for p,s in pssm.search(ts.seq,threshold=threshold_80perc,both=False): 
        if p >= flanking_size & p < Lpeak - Lmotif - flanking_size:
            # Write the site of the hit as the new locus
            seq = ts.seq[p - flanking_size : p+flanking_size+Lmotif]
            new_locus = p + peak_start
            newID = ':'.join((chr_name,str(new_locus),'+'))
            newrec = SeqRecord.SeqRecord( seq = seq, id = newID, description = ts.name )
            sites.append( newrec )
    
    rc = ts.reverse_complement()
    for p,s in pssm.search(rc.seq,threshold=threshold_80perc,both=False):
        if p >= flanking_size & p < Lpeak - Lmotif - flanking_size:
            # Write the "head" of the motif as the new locus
            seq = rc.seq[p - flanking_size : p+flanking_size+Lmotif]
            new_locus = peak_end - p
            newID = ':'.join((chr_name,str(new_locus),'-'))
            newrec = SeqRecord.SeqRecord( seq = seq, id = newID, description = ts.name )
            sites.append( newrec )

toc = default_timer()
SeqIO.write(sites,'/data/crispri_hamming/nanog_sites.fasta','fasta')

def get_peak_location(name):
    start = name.split(':')[1]
    end = start.split('-')[1]
    start = start.split('-')[0]
    return int(start), int(end)


def get_peak_chromosome(name):
    return name.split(':')[0]
