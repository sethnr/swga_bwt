#!/usr/bin/python

import sys
sys.path.remove('/usr/lib/python2.7/dist-packages')
sys.path.append('/usr/lib/python2.7/dist-packages')


from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil,factorial
import random as rand
import numpy as np
import h5py as h5
import os.path as path
from re import match
import profile
import itertools
import pickle as pkl

import swga_bwt as s


#############
# do some stuff
############

#backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
#           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']


if len(sys.argv)==1: 
  print "\n\tpython swga_bwt_gapfill.py index_file pattern_file no_patterns_to_check:index\n\n"
  sys.exit(1)

idxfile = sys.argv[1]
patternfile = sys.argv[2]
tmatchidx = sys.argv[3]
out = sys.stdout

pblocksize=-1
if len(sys.argv) > 4:
  pblocksize = sys.argv[4]
  pi = sys.argv[5]
  pblocksize = int(pblocksize)
  pi = int(pi)

#  out = open(outfile,'w')

pfile = open(patternfile,'r')
allPatterns = []
ratios = []
for line in pfile:
  if match('#',line):
    next
  else:
    F = line.split()
    if len(F) > 1:
      allPatterns += [F[0]]
      ratios += [F[-1]]
    else:
      allPatterns += [line.rstrip('\n').lower()]
      ratios += [1]

print allPatterns,"!"

# limit patterns to pi'th set of pblocksize
if pblocksize != -1:
  pstart = pblocksize * pi
  pend = pstart + pblocksize
  if pstart >= len(allPatterns): sys.exit(1)
  elif pend >= len(allPatterns): pend = len(allPatterns)-1 
  patterns = allPatterns[pstart:pend]

print >>sys.stderr, "indexing ", len(patterns), " of ",len(allPatterns)

      
#    print patterns
ratios = np.array(ratios,dtype='float')

#guess chrname and blocksize from indexname
fasta, blocksize = path.basename(idxfile).split('.')

blocksize = int(blocksize)
index = h5.File(idxfile+".IDX.hdf5", "r")
chrs = index.keys()

#initialise index if it doesn't exist
if not path.exists(tmatchidx):
#  tmatchfile = open(tmatchidx,"w")
  matchindex = h5.File(tmatchidx, "w")
#compression = None
  compression = 'gzip'
  matchindex =  s.initMatchCountsHDF5(index,allPatterns,blocksize,chrs,ratios,matchindex)
  print matchindex
  matchindex.close()


#update index with new patterns
#tmatchfile = open(tmatchidx,"r")
matchindex = h5.File(tmatchidx, "a")
s.addMatchCountsHDF5(index, patterns, matchindex)
matchindex.close()


