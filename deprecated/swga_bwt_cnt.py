#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
from re import match
#import profile


allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix


def firstColNP(tots):
  # get ranges of each base in first col (contiguous as BWM is sorted)
  first = {}
  totc = 0
  for n in range(0,len(allBases)):
    c = allBases[n]
    count = tots[n]
    first[c] = (totc, totc + count)
    totc += count
  return first

def countMatchesNP(baseRanks, first, p):
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

backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']

idxfile = sys.argv[1]
patternfile = sys.argv[2]
out_prefix = sys.argv[3]

#out = sys.stdout
#put proper argparser into this script asap

wig = open(out_prefix+".wig",'w')
out = open(out_prefix+".tab.txt",'w')

#if len(sys.argv) > 3:
#  outfile = sys.argv[3]
#  out = open(outfile,'w')


pfile = open(patternfile,'r')
patterns = []
patternset = []
pi = -1
for line in pfile:
  if match('#',line):
    next
  else:
    F = line.split()
    if len(F) > 1:
      if not match('[actgACTG,]',F[0]):        
        pi += 1
        allPatterns = F[-1]
        patterns = allPatterns.split(',')
        print >> sys.stderr, pi
        patternset.append( patterns)
#        print patternset
#        break
      else:
        patterns += [F[0]]
    else:
      patterns += [line.rstrip('\n').lower()]
if pi == -1:
  patternset[0] = [patterns]
#print patternset

#sys.exit(1)

#patterns = ['nnnnn'] + patterns
#guess chrname and blocksize from indexname
fasta, blocksize = path.basename(idxfile).split('.')

blocksize = int(blocksize)
#chr_index = np.load(idxfile+".IDX.npy")
#chr_bwts = np.load(idxfile+".BWT.npy")

index = h5.File(idxfile+".IDX.hdf5", "r")

#if len(sys.argv) > 3:
#  chrs = sys.argv[3:]
#else:
chrs = index.keys()



#for b in backgrounds:
#  background = h5file(b,"r")
#  b_index = index["subset/idx"]
#  b_bwts = index["subset/bwt"]
#  for p in patterns:
#    

pset = 0
for patterns in patternset:
  pset += 1
  wig = open(out_prefix+"."+str(pset)+".wig",'w')
  out = open(out_prefix+"."+str(pset)+".tab.txt",'w')

  uncovered = 0
  allBlocks = 0
  ptottot = 0
  for chr in chrs:
    chr_index = index[chr+"/idx"]
    chr_bwts = index[chr+"/bwt"]

    print >> sys.stderr, fasta, blocksize, chr
    blocks = chr_index.shape[0]
    
#sys.exit(1)

    print >> wig, "fixedStep chrom="+chr+" start=1 step="+str(blocksize)+" span="+str(blocksize)

    for n in range(0,blocks):
      allBlocks += 1 
      baseRanks = chr_index[n]
      bwt_line = chr_bwts[n]
  
  #get mapping of each base to range of posns in first column 
  #(all contiguous as col is sorted)
      firstColMap = firstColNP(baseRanks[-1])
      print >> out, chr, n*blocksize, (n+1)*blocksize, blocksize,
      ptot = 0
      for p in patterns:
        pCount =  countMatchesNP(baseRanks,firstColMap,p)
        ptot += pCount
        print >> out, pCount,
      print >> out, ""
    
      print >> wig, str(ptot)
      if(ptot == 0): uncovered += 1
      ptottot += ptot
  wig.close()
  out.close()
  print >>sys.stderr, "total blocks ("+str(blocksize)+"):"+str(allBlocks)
  print >>sys.stderr, "mean cov:"+str((ptottot / allBlocks))
  print >>sys.stderr, "unfilled:"+str(uncovered)
