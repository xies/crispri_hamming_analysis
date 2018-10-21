#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 15:19:45 2018

@author: xies
"""

import numpy as np
import re
import pandas as pd
from Bio import SeqIO, SeqRecord
from itertools import count, izip

# Load ChIP peaks
filename = '/data/crispri_hamming/oct4_chip_flanking.fa'
seqs = [rec for rec in SeqIO.parse(filename,'fasta')]
Npeaks = len(seqs)

# Load output of Ali's GG-finding file
guide_seqs = []
repeats = []
for rec in seqs:
    # Use positive lookahead to search through all matches
    pam_sites = [m.start() for m in re.finditer('(?=GG)',str(rec.seq))]
    # Filter through ones that have 20 flanking positions
    pam_sites = [s for s in pam_sites if s > 19]
    for s in pam_sites:
        peak_locus = get_genomic_locus_begin(rec.name)
        chr_name = get_chromosome_name(rec.name)
        new_locus = peak_locus + s - 21
        newID = ':'.join((chr_name,str(new_locus)))
        sequence = rec.seq[s-21:s-1]
        if len(sequence) == 0:
            break
        # Check for unique entries
        if newID not in [guide.id for guide in guide_seqs]:
            # Generate a SeqRecord
            guide_seqs.append( SeqRecord.SeqRecord(
                    seq = sequence,
                    id = newID,
                    name = rec.name) )
        else:
            repeats.append( SeqRecord.SeqRecord(
                    seq = rec.seq[s-21:s-1],
                    id = newID,
                    name = rec.name) )

SeqIO.write(guide_seqs,'/data/crispri_hamming/sgRNA_unique.fasta','fasta')
