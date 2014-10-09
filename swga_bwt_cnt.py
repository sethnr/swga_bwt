#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
#import profile


allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix


def firstColNP(tots):
    """ Return a map from characters to the range of cells in the first
column containing the character. """
    first = {}
    totc = 0
#    for c, count in sorted(tots.iteritems()):
    for n in range(0,len(allBases)):
        c = allBases[n]
        count = tots[n]
        first[c] = (totc, totc + count)
        #first[n] = (totc, totc + count)
        totc += count
    return first

def countMatchesNP(baseRanks, first, p):
    """ Given BWT(T) and a pattern string p, return the number of times
p occurs in T. """
    if p[-1] not in first:
        return 0 # character doesn't occur in T
    l, r = first[p[-1]]
    i = len(p)-2
    while i >= 0 and r > l:
        c = p[i]
        ci = allBases.index(c)
        l = first[c][0] + baseRanks[l-1,ci]
        r = first[c][0] + baseRanks[r-1,ci]
        i -= 1
    return r - l # return size of final range


#############
# do some stuff
############

blocksize = 5000

idxfile = sys.argv[1]
patternfile = sys.argv[2]
pfile = open(patternfile,'r')
patterns = []
for line in pfile:
  F = line.split()
  patterns += [F[0]]
#print patterns

#sys.exit(1)
fasta, chrname, blocksize = path.basename(idxfile).split('.')
blocksize = int(blocksize)

chr_index = np.load(idxfile+".IDX.npy")
chr_bwts = np.load(idxfile+".BWT.npy")

print fasta, blocksize, chrname
blocks = chr_index.shape[0]

#sys.exit(1)

for n in range(0,blocks):
  baseRanks = chr_index[n]
  bwt_line = chr_bwts[n]
  
  #get mapping of each base to range of posns in first column 
  #(all contiguous as col is sorted)
  firstColMap = firstColNP(baseRanks[-1])
  print n, n*blocksize, (n+1)*blocksize, blocksize,
  for p in patterns:
    print countMatchesNP(baseRanks,firstColMap,p),
  print ""
