#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 15:19:45 2018

@author: xies
"""

import numpy as np
import re
from Bio import SeqIO, SeqRecord, Seq
from itertools import product

# Load ChIP peaks
filename = '/data/crispri_hamming/oct4_chip_flanking.fa'
peaks = [rec for rec in SeqIO.parse(filename,'fasta')]
Npeaks = len(peaks)
has_guide = np.zeros(Npeaks,dtype=bool)
Lpeak = len(peaks[0]) - 1

# Load output of Ali's GG-finding file
guide_seqs = []
repeats = []
for (i,rec) in enumerate(peaks):
    # Use positive lookahead to search through all matches on + strand
    pam_sites_plus = list(product([m.start() for m in re.finditer('(?=GG)',str(rec.seq))],'+'))
    # Combine with the reverse complement to search through on - strand
    revcomp = str(rec.reverse_complement().seq)
    pam_sites_minus = list(product([m.start() for m in re.finditer('(?=GG)',revcomp)],'-'))
    # Filter through ones that have 20 flanking positions
    pam_sites = [s for s in pam_sites_plus if s[0] > 19]
    pam_sites = pam_sites + [s for s in pam_sites_minus if s[0] > 19]
    for s in pam_sites:
        strand = s[1]; s = s[0]
        chr_name = get_peak_chromosome(rec.name)
        # Grab sequence depending on strand
        if strand == '+':
            peak_locus = get_peak_location(rec.name)[0]
            new_locus = peak_locus + s - 21
            sequence = rec.seq[s-21:s-1]
            newID = ':'.join((chr_name,str(new_locus),'+'))
        else:
            peak_locus = get_peak_location(rec.name)[1]
            new_locus = peak_locus - s + 21
            newID = ':'.join((chr_name,str(new_locus),'-'))
            sequence = revcomp[s-21:s-1]
            
        if len(sequence) == 0:
            break
        # Check for unique entries
        if newID not in [guide.id for guide in guide_seqs]:
            if strand == '+':
                sequence = rec.seq[s-21:s-1]
                newrec = SeqRecord.SeqRecord( seq = sequence, id = newID, description = rec.name )
            elif strand == '-':
                sequence = Seq.Seq( data=revcomp[s-21:s-1] )
                newrec_rc = SeqRecord.SeqRecord( seq = sequence )
                newrec.id = newID
                newrec.description = rec.name
            # Generate a SeqRecord
            guide_seqs.append( newrec)
            # note that this peak has a guide
            has_guide[i] = True
        else:
            has_guide[i] = True
            repeats.append( SeqRecord.SeqRecord(
                    seq = rec.seq[s-21:s-1],
                    id = newID,
                    description = rec.name) )
    print "Done with ", i

print float(has_guide.sum()) / Npeaks

SeqIO.write(guide_seqs,'/data/crispri_hamming/sgRNA_unique.fasta','fasta')
#0.90838375108

def get_peak_location(name):
    start = name.split(':')[1]
    end = start.split('-')[1]
    start = start.split('-')[0]
    return int(start), int(end)

def get_peak_chromosome(name):
    return name.split(':')[0]

