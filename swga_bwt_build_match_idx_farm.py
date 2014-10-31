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
  print "\n\tpython swga_bwt_build_match_idx.py index_file pattern_file no_patterns_to_check:index\n\n"
  sys.exit(1)

#parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
#parser.add_argument('-t|--target', dest='idxfile', type=str, help='target genome index')
#parser.add_argument('-p|--patterns', dest='patternfile', type=strhelp='list of patterns')
#parser.add_argument('-p|--tmatchidx', dest='tmatchidx', type=str, help='match positions index')
#parser.add_argument('-o|--out', dest='out', type=str, help='outfile')

#args = parser.parse_args()

pi = os.getenv("LSB_JOBINDEX")

idxfile = sys.argv[1]
patternfile = sys.argv[2]
patternfile = patternfile+"."+pi

tmatchidx = sys.argv[3]
out = sys.stdout


#  out = open(outfile,'w')

pfile = open(patternfile,'r')
patterns = []
ratios = []
for line in pfile:
  if match('#',line):
    next
  else:
    F = line.split()
    if len(F) > 1:
      patterns += [F[0]]
      ratios += [F[-1]]
    else:
      patterns += [line.rstrip('\n').lower()]
      ratios += [1]

      
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


