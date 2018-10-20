#!/usr/bin/python


"""
Finds the pairwise Hamming distance between each FASTA sequence
in the input file.
"""

import numpy as np
import scipy as sp
from Bio import SeqIO, AlignIO
import sys, getopt, timeit
from itertools import izip, count, product 

def main(argv):
	# Parse some inputs
#	try:
#		opts,args = getopt.getopt(argv,'hi:o',['ifile=','ofile='])
#	except getopt.GetoptError:
#		print 'pairwise_hamming_distance.py -i <inputfile> -o <outputfile>'
#		sys.exit(2)
#	for opt,args in opts:
#		if opt == '-h':
#			print 'pairwise_hamming_distance.py -i <inputfile> -o <outputfile>'
#			sys.exit(2)
#		elif opt in ('-i','--ifile'):
#			inputfile = arg
#		elif opt in ('-o','--ofile'):
#			outputfile = arg

	inputfile = argv[0]

	# Load file
	seqs = list(SeqIO.parse(inputfile,'fasta'))
	seqs = [str(x.seq) for x in seqs]
	Ntotal_seqs = len(seqs)

	# Get the subset of sequences
	iterbegin = int(argv[1])
	iterend = int(argv[2])

	# Check that the subset are within length
	if iterbegin + iterend - 1 < Ntotal_seqs:
		seqSub = seqs[iterbegin - 1 : (iterbegin + iterend - 1)]
	else:
		seqSub = seqs[iterbegin - 1:]

	# Build the index of pairs
	I = xrange(len(seqSub))
	J = xrange(Ntotal_seqs) # needs to be changed for subsetting
	index = list(product(I,J))
		
	D = np.zeros((len(seqSub),len(seqs)),dtype=np.int)
	# Iterate through pariwise lists using itertools
	for i, seqPair in izip(count(), product(seqSub,seqs)):
		ind = index[i]
		# If along the diagonal, put 100 instead (to not interfere with min)
		if (ind[0] + iterbegin-1) == ind[1]:
			D[ind[0],ind[1]] = 100
		else:
			D[ind[0],ind[1]] = hamming_distance(seqPair[0],seqPair[1])
	
	# Take the minimum distance except along diagonal (which will never be min)
	Dargmin = D.argmin(axis=1)
	D = D.min(axis=1)
	D = np.vstack( (D,Dargmin) ).T

	outputfile = argv[3]
	np.savetxt(outputfile, D.astype(np.int), fmt='%i', delimiter=',')

def hamming_distance(s1,s2):
	return sum(x != y for x,y in zip(s1,s2))

if __name__ == '__main__':
	main(sys.argv[1:])
	
