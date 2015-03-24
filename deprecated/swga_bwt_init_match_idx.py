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
import os
from re import match
import profile
import itertools
import pickle as pkl 
import argparse
import swga_bwt as s


#############
# do some stuff
############

#backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
#           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']


if len(sys.argv)==1: 
  print "\n\tpython swga_bwt_build_match_idx.py index_file pattern_file no_patterns_to_check:index\n\n"
  sys.exit(1)

parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-t|--target', action="store", dest='idxfile', type=str, help='target genome index', nargs='?')
parser.add_argument('-p|--patterns', action="store", dest='patternfile', type=str, help='list of patterns', nargs='?')
parser.add_argument('-i|--tmatchidx', action="store", dest='tmatchidx', type=str, help='match positions index', nargs='?')
parser.add_argument('-o|--out', action="store", dest='out', type=str, help='outfile', nargs='?')

args = parser.parse_args()

#pi = os.getenv("LSB_JOBINDEX")

print 
idxfile = args.idxfile
patternfile = args.patternfile
#patternfile = patternfile+"."+pi

tmatchidx = args.tmatchidx
out = args.out
if out == None:
  out = sys.stdout
else:
  out = open(outfile,'w')



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
  matchindex = h5.File(tmatchidx, "w")
  compression = 'gzip'
  matchindex =  s.initMatchCountsHDF5(index,patterns,blocksize,chrs,ratios,matchindex)
  print matchindex
  matchindex.close()
