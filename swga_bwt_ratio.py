#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
#import profile
import pickle as pkl
from swga_bwt import *


#############
# do some stuff
############
target = '/nfs/users/nfs_s/snr/swga/idx/Pf3D7_v3.0.1.3000'
backgrounds = {'moz':'/nfs/users/nfs_s/snr/swga/idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.0.01.3000',
               'man':'/nfs/users/nfs_s/snr/swga/idx/Homo_sapiens.GRCh38.dna_sm.primary_assembly.0.001.3000'}


if len(sys.argv)==1: 
  print "\n\tpython swga_bwt_ratio.py pattern_file outfile countindex no_patterns_to_check:index\n\n"
  sys.exit(1)

#idxfile = sys.argv[1]
patternfile = sys.argv[1]

#out = sys.stdout
if len(sys.argv) > 2:
  outfile = sys.argv[2]
  if outfile == '-':
    out = sys.stdout
  else:
    out = open(outfile,'w')

if len(sys.argv) > 3:
  idxfile = sys.argv[3]

pblocksize = -1
if len(sys.argv) > 4:
  pblocksize,pi = sys.argv[4].split(':')
pblocksize = int(pblocksize)
pi = int(pi)

print >>sys.stderr, patternfile
print >>sys.stderr, outfile
print >>sys.stderr, idxfile
print >>sys.stderr, "noPats", pblocksize, "index",pi


#guess chrname and blocksize from indexname
filebits = path.basename(target).split('.')
fasta, blocksize = filebits[0],filebits[-1]
blocksize = int(blocksize)
blockToRPK = float(1000) / blocksize
#print >> sys.stderr, "block to RPK:",blockToRPK

pfile = open(patternfile,'r')
patterns = []
tcounts={}

#read in pattern file
for line in pfile:
  pattern, count, countBlock = line.split()
  #print "pattern",pattern,"count",count,"RPK",countBlock
  patterns += [pattern]
  tcounts[pattern] = float(countBlock) #* blockToRPK

allPatterns = patterns #save this, we'll need it later

# limit patterns to pi'th set of pblocksize
if pblocksize != -1:
  pstart = pblocksize * pi
  pend = pstart + pblocksize
  if pstart >= len(patterns): sys.exit(1)
  elif pend >= len(patterns): pend = len(patterns)-1 
  patterns = patterns[pstart:pend]

print >>sys.stderr, len(patterns)


#blocksize = int(blocksize)
#chr_index = np.load(idxfile+".IDX.npy")
#chr_bwts = np.load(idxfile+".BWT.npy")

#index = h5.File(target+".IDX.hdf5", "r")


if path.exists(idxfile):
  cachecounts = open(idxfile,"r")
  backcounts = pkl.load(cachecounts)
  cachecounts.close()

  indexPatterns = [str(i[0]) for i in backcounts.keys()]
  newPatterns = [p for p in patterns if p not in indexPatterns]
  print >> sys.stderr, "patterns:", len(patterns), len(newPatterns)
  if(len(newPatterns) > 0):
    newCounts = getIndexCounts(backgrounds,newPatterns, blockToRPK)
    backcounts.update(newCounts)
    cachecounts = open(idxfile, 'w')
    pkl.dump(backcounts,cachecounts)
    cachecounts.close()
else:
  backcounts = getIndexCounts(backgrounds, patterns, blockToRPK)
  cachecounts = open(idxfile,"w")
  pkl.dump(backcounts,cachecounts)
  cachecounts.close()
#print backcounts

print >> out, "#primer", "t_count",
for b in backgrounds.keys():
  print >>out, str(b)+"_count", str(b)+"_ratio",
print >> out, "mean_ratio"

for p in patterns:
#  print p, target
#  tcount = backcounts[p,target+".IDX.hdf5"]
  tcount = tcounts[p]
  ratioTotal = 0
  print >>out, p, tcount, 
  for b in backgrounds.values():
    bcount = backcounts[p,b+".IDX.hdf5"]
    if bcount > 0: ratio = tcount/bcount
    else: ratio=float('inf')
    print >>out, bcount, ratio,
    ratioTotal += ratio
  print >>out, ratioTotal / len(backgrounds)
#  print >>out, ''

sys.exit(1)
