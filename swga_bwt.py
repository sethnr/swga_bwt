import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
import pickle as pkl
from re import findall

allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix


def getTM(primer):
    nA = len(findall('a',primer))
    nC = len(findall('c',primer))
    nT = len(findall('t',primer))
    nG = len(findall('g',primer))
    tm = nA*2+nT*2+nC*4+nG*4
    #print tm
    return tm




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
#  for b in backgrounds.values():
  for b in backgrounds:
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
        backcounts[(p,b)] = (float(matches) / float(noBlocks)) * blockToRPK
      print '.',
      #print p, matches, noBlocks
  return backcounts


def initMatchCountsHDF5(index, allPatterns, blocksize, chrs, ratios,matchindex):
  allBlocks = 0
  for chr in chrs:
    chr_bwts = index[chr+"/bwt"]
    allBlocks += chr_bwts.shape[0]
  matchcountsNP = np.zeros(shape=(allBlocks,len(allPatterns)),dtype='int')
  blockPosnsNP = np.zeros(shape=(allBlocks,2),dtype='int')
  blockChrsNP = np.empty(shape=(allBlocks,1),dtype='S25')
  # print >>sys.stderr, allPatterns
 
  bi = -1
  for chr in chrs:
    chr_index = index[chr+"/idx"]
    blocks = chr_index.shape[0]
    for n in range(0,blocks):
      bi += 1
      blockChrsNP[bi] = chr
      blockPosnsNP[bi] = (n*blocksize,(n+1)*blocksize)

#  pArray = np.recarray(shape=(len(allPatterns),2),dtype=[('pattern', 'S20'), ('i', int)])
#  print >>sys.stderr, allPatterns
#  pArray['pattern'] = allPatterns
#  pArray['i'] = range(0,len(allPatterns))
  
  
#  matchindex['matchcounts'] = matchcountsNP
#  matchindex['blockposns'] = blockPosnsNP
#  matchindex['patterns'] = allPatterns
#  matchindex['ratios'] = np.array(ratios,dtype='float')


  matchcounts = matchindex.create_dataset('matchcounts', 
                                          (allBlocks,len(allPatterns)), 
                                          dtype='int',
                                          data=matchcountsNP)
  blockposns = matchindex.create_dataset('blockposns', 
                                         (allBlocks,2), 
                                         dtype='int',
                                         data=blockPosnsNP)
  blockchrs = matchindex.create_dataset('blockchrs', 
                                        (allBlocks,), 
                                        dtype='S25',
                                        data=blockChrsNP)
  patterns = matchindex.create_dataset('patterns', 
                                       (len(allPatterns),), 
                                       dtype='S20',
                                       data = allPatterns)  
  ratios = matchindex.create_dataset('ratios', 
                                     (len(allPatterns),), 
                                     dtype='float',
                                     data=ratios)

  print patterns[:]
  return matchindex

#non numpy version, so we can use patterns as indices:
def getMatchCounts(index, patterns, blocksize, chrs):
  
#  allBlocks = matchcounts.shape[0]
  matchcounts = {}
  blockPosns = {}
#  matchcounts = np.zeros(shape=(allBlocks,len(patterns)),dtype='int')
#  blockPosns = np.recarray(shape=(allBlocks,3),dtype=[('chr', 'S12'), ('start', int),('end',int)])
  bi=-1
  for chr in chrs:
    chr_index = index[chr+"/idx"]
    chr_bwts = index[chr+"/bwt"]
    
    print >> sys.stderr, chr, chr_index.shape, chr_index, blocksize
    blocks = chr_index.shape[0]
    
#    print >>sys.stderr, "blocks = ",blocks
    for n in range(0,blocks):
      baseRanks = chr_index[n]
      bwt_line = chr_bwts[n]
      bi +=1
    #get mapping of each base to range of posns in first column 
    #(all contiguous as col is sorted)
      firstColMap = firstColNP(baseRanks[-1])
    
#      print >> sys.stderr, chr, n*blocksize, (n+1)*blocksize, blocksize, bi
      
      blockPosns[bi] = (chr,n*blocksize,(n+1)*blocksize)
      #blockPosns[(bi,0:2)] = (chr,n*blocksize,(n+1)*blocksize)
      for p in patterns:
        matchcounts[(bi,p)] = countMatchesNP(baseRanks,firstColMap,p)
  return matchcounts, blockPosns

#numpy version
def addMatchCountsHDF5(index, newPatterns, matchcountsHDF5):
  
  allBlocks = matchcountsHDF5['matchcounts'].shape[0]
  matchcounts = matchcountsHDF5['matchcounts']
  blockPosns = matchcountsHDF5['blockposns']
  blockChrs = matchcountsHDF5['blockchrs'][:]
  allPatterns = matchcountsHDF5['patterns'][:]
  
#  matchcounts = np.zeros(shape=(allBlocks,len(patterns)),dtype='int')
#  blockPosns = np.recarray(shape=(allBlocks,3),dtype=[('chr', 'S12'), ('start', int),('end',int)])

  # GET POSITIONS OF PATTERNS IN FILE
  # for p_index in P_indices
  
  patterns = allPatterns[np.in1d(allPatterns, newPatterns)]
  
#  pis = allPatterns[np.in1d(allPatterns.pattern, newPatterns),'i']
  
  print >>sys.stderr, patterns
#  print >>sys.stderr, pis
  
  chrs = set(blockChrs)
  bi=-1
  for chr in chrs:
    print >> sys.stderr, chr,
    chr_index = index[chr+"/idx"]
    chr_bwts = index[chr+"/bwt"]
    
    print >> sys.stderr, chr_index.shape, chr_index
    blocks = chr_index.shape[0]
    
#    print >>sys.stderr, "blocks = ",blocks
    for n in range(0,blocks):
      baseRanks = chr_index[n]
      bwt_line = chr_bwts[n]
      bi +=1
    #get mapping of each base to range of posns in first column 
    #(all contiguous as col is sorted)
      firstColMap = firstColNP(baseRanks[-1])
    
#      print >> sys.stderr, chr, n*blocksize, (n+1)*blocksize, blocksize, bi
      
      #blockPosns[bi] = (chr,n*blocksize,(n+1)*blocksize)
      #blockPosns[bi,0:2] = (chr,n*blocksize,(n+1)*blocksize)
      for p in patterns:
        pi = int(np.where(allPatterns==p)[0])
#        print >>sys.stderr, p, pi, bi
        matchcounts[bi,pi] = countMatchesNP(baseRanks,firstColMap,p)
  


#def updateMatchCounts(matchcounts, blockPosns, ):
#  patterns = list(matchcounts.keys(), key=lambda a: a[2])
#
#  which = np.in1d(pArray['pattern'], patterns)
#  pis = parray['i'][which]
#  patternsNP = parray['pattern'][which]
#  for i in range(0,len(pis)):
#    matchcountsNP[bi,pis[i]] = matchcounts[(bi,pis[i])]
#  print >>sys.stderr, p,"->",pi
  

#  return matchcounts, blockPosns #, pArray
#    print >> out, ""


########
# BWT index builder
########

def suffixArray(s):
  #get tuples of suffixes + ordering of original rows
  sa = sorted([(s[i:], i) for i in xrange(0, len(s)+1)])
  #return just ordering of rows
  return map(lambda x: x[1], sa)

def bwt(t):
  #get burrows wheeler transform of string T
  bw = []
  for si in suffixArray(t):
    if si == 0:
      bw.append('$')
    else:
      bw.append(t[si-1])
  return ''.join(bw) # return string-ized version of list bw


def rankAllBwtNP(bw):
  # returns an ndarray of size (seq_length * no_different_bases) 
  # containing cumulative nos of occurences of base n
  tots_i = np.array([0]*len(allBases))
  i = 0;
  baseRanks = np.zeros(shape=(len(bw),len(allBases)),dtype='int')

  for c in bw:
      c = c.lower()      
      if c not in allBases:
          c = 'n'
      nc = allBases.index(c)
      tots_i[nc] += 1
      baseRanks[i] = tots_i
      i +=1
  return baseRanks
