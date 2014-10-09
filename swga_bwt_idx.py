#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
#import profile


allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix

def suffixArray(s):
  """ Given T return suffix array SA(T). We use Python's sorted
function here for simplicity, but we can do better. """
  satups = sorted([(s[i:], i) for i in xrange(0, len(s)+1)])
  # Extract and return just the offsets
 # print satups
  return map(lambda x: x[1], satups)

def bwt(t):
  """ Given T, returns BWT(T), by way of the suffix array. """
  bw = []
  for si in suffixArray(t):
    if si == 0:
      bw.append('$')
    else:
      bw.append(t[si-1])
  return ''.join(bw) # return string-ized version of list bw


def rankAllBwtNP(bw):
  """ Given BWT string bw, returns a map of lists. Keys are
characters and lists are cumulative # of occurrences up to and
including the row. """

  tots_i = np.array([0]*len(allBases))
  #tots_r = np.rec.array([0]*len(allBases))
  #tots_r.dtype.names = allBases
  i = 0;
  baseRanks = np.zeros(shape=(len(bw),len(allBases)),dtype='int')
#  baseRanks = np.ndarray(shape=(len(bw),len(allBases)),dtype=int)
  for c in bw:
      c = c.lower()      
      if c not in allBases:
          c = 'n'
      nc = allBases.index(c)
      tots_i[nc] += 1
      baseRanks[i] = tots_i
      #tots_r[c] += 1
      #baseRanks[i] = tots_r
      i +=1
  return baseRanks

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

def getTotsNP(baseRanks):
  print "WTF!"


def countMatchesNP(baseRanks, first, p):
    """ Given BWT(T) and a pattern string p, return the number of times
p occurs in T. """
#    baseRanks = rankAllBwtNP(bw)
    rankAll = {}
#    tots = {}
#    for i in range(0,len(allBases)):
#        rankAll[allBases[i]] = baseRanks[:,i].tolist()
#        tots[allBases[i]] = baseRanks[-1,i]

#    tots = baseRanks[-1]
#    first = firstColNP(tots)
    if p[-1] not in first:
        return 0 # character doesn't occur in T
    l, r = first[p[-1]]
    i = len(p)-2
    while i >= 0 and r > l:
        c = p[i]
        ci = allBases.index(c)
        #l = first[ci][0] + baseRanks[l-1,ci]
        #r = first[ci][0] + baseRanks[r-1,ci]
        l = first[c][0] + baseRanks[l-1,ci]
        r = first[c][0] + baseRanks[r-1,ci]
        i -= 1
    return r - l # return size of final range


#############
# do some stuff
############

#blocksize = 5000

fasta = sys.argv[1]
blocksize = int(sys.argv[2])

print fasta
seq = SeqIO.parse(open(fasta),'fasta')

p = 'act'

ci = 0
for chr in seq:
  ci += 1
  if ci > 1: sys.exit(1)
  name = chr.id    
  seqlen = len(chr.seq)
  print name, seqlen
  lastBlock = (seqlen/blocksize)
    #lastBlock #rounds down as both are integers
  
  #indices - add one to blocksize for EOL marker ($)
  chr_index = np.empty(shape=(lastBlock+1,blocksize+1,len(allBases)),dtype='int')
  chr_bwts = np.empty(shape=(lastBlock+1,blocksize+1),dtype='string_')

  for n in range(0,lastBlock+1):
    end = (n+1)*blocksize
    if end > seqlen: #if not full block, pad with Ns
      pad = 'n'*(end-seqlen)
#      print len(pad)
      seqblock = chr.seq[n*blocksize:seqlen]
      seqblock = seqblock + Seq(pad)
    else:  
      seqblock = chr.seq[n*blocksize:end]
      
    print n, n*blocksize, end, len(seqblock)
    bwt_line = bwt(str(seqblock))     
    baseRanks = rankAllBwtNP(bwt_line)
    
    #get mapping of each base to range of posns in first column 
    #(all contiguous as col is sorted)
#    firstColMap = firstColNP(baseRanks[-1])
#    print firstColMap
#    print seqblock.count(p),
#    print countMatchesNP(baseRanks,firstColMap,p)
     
    chr_index[n] = baseRanks
    chr_bwts[n] = bwt_line
#    print countMatchesNP(baseRanks,p)        
#        profile.run("rankAllBwtNP(bwt_line)")
  fasta = fasta.replace('.fasta','')
  fasta = fasta.replace('.fa','')
  idxfile = "./idx/"+fasta+"."+name+"."+str(blocksize)
#  idxfile = fasta+"."+name+"."+str(blocksize)
  np.save(idxfile+".IDX.npy",chr_index)
  np.save(idxfile+".BWT.npy",chr_bwts)  
    #if n > 10: sys.exit(1)

