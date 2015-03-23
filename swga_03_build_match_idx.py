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
import argparse
import os

#############
# do some stuff
############


if len(sys.argv)==1: 
  print "\n\tpython swga_bwt_build_match_idx.py -t target_index -p pattern_file -i match_index \n\n"
  sys.exit(1)


parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-t','--target', action="store", dest='idxfile', type=str, help='target genome index', nargs='?')
parser.add_argument('-p','--patterns', action="store", dest='patternfile', type=str, help='list of patterns', nargs='?')
parser.add_argument('-i','--tmatchidx', action="store", dest='tmatchidx', type=str, help='match positions index', nargs='?')
parser.add_argument('-o','--out', action="store", dest='out', type=str, help='outfile', nargs='?')
parser.add_argument('-A','--array', action="store_true", dest='farm', help='is this being run on the farm? i.e. use LSB_JOBINDEX')
parser.add_argument('-I','--init', action="store_true", dest='init', help='initialise matrix (will be done automatically if local)')

args = parser.parse_args()

#pi = os.getenv("LSB_JOBINDEX")

farm = args.farm
idxfile = args.idxfile
idxfile = idxfile.replace(".IDX.hdf5","")
patternfile = args.patternfile
tmatchidx = args.tmatchidx

#print >>sys.stderr, args.farm

if args.farm is True:
  pi = os.getenv("LSB_JOBINDEX")
  if pi == None:
    sys.exit(100,"LSB_JOBINDEX = none - did you mean to submit an array job?")
  patternfile = patternfile+"."+pi

out = args.out
if out == None:
  out = sys.stdout
else:
  out = open(outfile,'w')

  

pfile = open(patternfile,'r')
patterns = []
ratios = []
for line in pfile:
  if match('#',line): next
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


def _initialiseMatchIDX():
  matchindex = h5.File(tmatchidx, "w")
  compression = 'gzip'
  matchindex =  s.initMatchCountsHDF5(index,patterns,blocksize,chrs,ratios,matchindex)
  print matchindex
  matchindex.close()

## initialise match array
#NB potential race condition with matrix initialisation in farm use
# therefore initialise only if local job or specific init job
if not path.exists(tmatchidx):
    if args.farm is False:
        _initialiseMatchIDX()
    elif args.init is True:
        _initialiseMatchIDX()
    
#update index with new patterns
#tmatchfile = open(tmatchidx,"r")
# open match file for appending
matchindex = h5.File(tmatchidx, 'r+')
s.addMatchCountsHDF5(index, patterns, matchindex)
matchindex.close()


