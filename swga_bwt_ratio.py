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
  return int(r - l) # return size of final range


#############
# do some stuff
############
target = './idx/Pf3D7_v3.0.1.3000'
backgrounds = {'moz':'./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.0.01.3000',
               'man':'./idx/Homo_sapiens.GRCh38.dna_sm.primary_assembly.0.001.3000'}

#idxfile = sys.argv[1]
patternfile = sys.argv[1]

out = sys.stdout

if len(sys.argv) > 2:
  outfile = sys.argv[2]
  out = open(outfile,'w')

#guess chrname and blocksize from indexname
filebits = path.basename(target).split('.')
fasta, blocksize = filebits[0],filebits[-1]
blocksize = int(blocksize)
blockToRPK = float(1000) / blocksize
print >> sys.stderr, "block to RPK:",blockToRPK

pfile = open(patternfile,'r')
patterns = []
tcounts={}
for line in pfile:
  #F = line.split()
  #patterns += [F[0]]
  pattern, count, countBlock = line.split()
  print "pattern",pattern,"count",count,"RPK",countBlock
  patterns += [pattern]
  tcounts[pattern] = float(countBlock) #* blockToRPK


#blocksize = int(blocksize)
#chr_index = np.load(idxfile+".IDX.npy")
#chr_bwts = np.load(idxfile+".BWT.npy")

index = h5.File(target+".IDX.hdf5", "r")

if len(sys.argv) > 3:
  chrs = sys.argv[3:]
else:
  chrs = index.keys()

def getIndexCounts(target, backgrounds, patterns):
  backcounts = {}
  print backgrounds
  #get target
  #tindex = h5.File(str(target)+".IDX.hdf5","r")
  #t_index = tindex["subset/idx"]
  #t_bwts = tindex["subset/bwt"]
  
  #for b in [target] + backgrounds.values():
  for b in backgrounds.values():
    print b
    b = str(b)+".IDX.hdf5"
    print b
    bgindex = h5.File(b,"r")
    b_index = bgindex["subset/idx"]
    b_bwts = bgindex["subset/bwt"]
  #firstColMap = firstColNP(baseRanks[-1])
    i = 0
    noBlocks = b_bwts.shape[0]
    print >> sys.stderr, b, noBlocks
    
    for p in patterns:
      i += 1
      if i % 10 == 0: print >>sys.stderr, i 
    #print >> out, chr, n*blocksize, (n+1)*blocksize, blocksize,
#      noBlocks = b_bwts.shape[0]
      #noBlocks = 10
      matches = 0
#      print p,
      for n in range(0,noBlocks):
        baseRanks = b_index[n]
        firstColMap = firstColNP(baseRanks[-1])
        bwt_line = b_bwts[n]
        match = countMatchesNP(baseRanks,firstColMap,p)
#        print match, noBlocks,
        matches += match
        backcounts[p,b] = (float(matches) / float(noBlocks)) * blockToRPK
      print p, matches, noBlocks
  return backcounts

if path.exists("index_counts.pkl"):
  cachecounts = open("index_counts.pkl","r")
  backcounts = pkl.load(cachecounts)
  indexPatterns = [str(i[0]) for i in backcounts.keys()]
  newPatterns = [p for p in patterns if p not in indexPatterns]
  print "patterns:", len(patterns), len(newPatterns)
  if(len(newPatterns) > 0):
    cachecounts.close()
    cachecounts = open("index_counts.pkl", 'w')
    newCounts = getIndexCounts(target,backgrounds,newPatterns)
    backcounts.update(newCounts)
    pkl.dump(backcounts,cachecounts)
    cachecounts.close()
else:
  backcounts = getIndexCounts(target, backgrounds, patterns)
  cachecounts = open("index_counts.pkl","w")
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
