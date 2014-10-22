import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
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

# get counts of no of times each pattern is represented in each background
def getIndexCounts(backgrounds, patterns, blockToRPK):
  backcounts = {}
  print backgrounds
  #get target
  #tindex = h5.File(str(target)+".IDX.hdf5","r")
  #t_index = tindex["subset/idx"]
  #t_bwts = tindex["subset/bwt"]
  
  #for b in [target] + backgrounds.values():
  for b in backgrounds.values():
    #print b
    b = str(b)+".IDX.hdf5"
    #print b
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
      print '.',
      #print p, matches, noBlocks
  return backcounts


def initiateMatchCounts(index, allPatterns, chrs):
  allBlocks = 0
  for chr in chrs:
    chr_bwts = index[chr+"/bwt"]
    allBlocks += chr_bwts.shape[0]
  matchcounts = np.zeros(shape=(allBlocks,len(allPatterns)),dtype='int')
  blockPosns = np.recarray(shape=(allBlocks,3),dtype=[('chr', 'S12'), ('start', int),('end',int)])
  blockPosns.fill('0')
  pArray = np.recarray(shape=(len(allPatterns),2),dtype=[('pattern', 'S20'), ('i', int)])
  parray['pattern'] = allPatterns
  parray['i'] = range(0,len(allPatterns))
  return matchcounts, blockPosns, pArray


def getMatchCounts(index, matchcounts, blockposns, pArray, patterns, blocksize, chrs):
  
  allBlocks = matchcounts.shape[0]

#  matchcounts = np.zeros(shape=(allBlocks,len(patterns)),dtype='int')
#  blockPosns = np.recarray(shape=(allBlocks,3),dtype=[('chr', 'S12'), ('start', int),('end',int)])

#def getMatchCounts(index, chrs):
  #sys.exit(1)

  which = np.in1d(pArray, patterns)
  pis = parray['i'][which]
  allPatterns = parray['pattern'][which]

  bi=-1
  for chr in chrs:
    chr_index = index[chr+"/idx"]
    chr_bwts = index[chr+"/bwt"]

    print >> sys.stderr, fasta, blocksize
    blocks = chr_index.shape[0]
    
    for n in range(0,blocks):
      baseRanks = chr_index[n]
      bwt_line = chr_bwts[n]
      bi +=1
    #get mapping of each base to range of posns in first column 
    #(all contiguous as col is sorted)
      firstColMap = firstColNP(baseRanks[-1])
    
    #    print >> out, chr, n*blocksize, (n+1)*blocksize, blocksize,
      
      blockPosns[bi,0:2] = (chr,n*blocksize,(n+1)*blocksize)
      for pi in pis:
        p = allPatterns[pi]
        matchcounts[bi,pi] = countMatchesNP(baseRanks,firstColMap,p)
#      print >> out, countMatchesNP(baseRanks,firstColMap,p),
  return matchcounts, blockPosns, pArray
#    print >> out, ""


