#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
#import profile


allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix

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

#############
# do some stuff
############

#blocksize = 5000

fasta = sys.argv[1]
blocksize = int(sys.argv[2])

print fasta
seq = SeqIO.parse(open(fasta),'fasta')

fasta = path.basename(fasta)
fasta = fasta.replace('.fasta','')
fasta = fasta.replace('.fa','')
idxfile = "./idx/"+fasta+"."+str(blocksize)

index = h5.File(idxfile+".IDX.hdf5", "w")
#compression = None
compression = 'gzip'

#ci = 0
for chr in seq:
#  ci += 1
#  if ci > 1: sys.exit(1)

  
  name = chr.id    
  seqlen = len(chr.seq)
#  print name, seqlen, blocksize, lastBlock
  lastBlock = (seqlen/blocksize)
    #lastBlock #rounds down as both are integers
  print name, seqlen, blocksize, lastBlock
  
  #indices - add one to blocksize for EOL marker ($)
#  chr_index = np.empty(shape=(lastBlock+1,blocksize+1,len(allBases)),dtype='int')
#  chr_bwts = np.empty(shape=(lastBlock+1,blocksize+1),dtype='string_')

  chr_index = index.create_dataset(name+"/idx", 
                                   (lastBlock+1,blocksize+1,len(allBases)), 
                                   dtype='i',
                                   compression=compression)
  chr_bwts = index.create_dataset(name+"/bwt", 
                                  (lastBlock+1,blocksize+1), 
                                  dtype='S1',
                                  compression = compression)
#  print chr_bwts.compression


  for n in range(0,lastBlock+1):
    end = (n+1)*blocksize
    if end > seqlen: #if not full block, pad with Ns
      pad = 'n'*(end-seqlen)
      seqblock = chr.seq[n*blocksize:seqlen]
      seqblock = seqblock + Seq(pad)
    else:  
      seqblock = chr.seq[n*blocksize:end]
      
#    print n, n*blocksize, end, len(seqblock)
    bwt_line = bwt(str(seqblock))     
    baseRanks = rankAllBwtNP(bwt_line)
    
    chr_index[n] = baseRanks
    chr_bwts[n] = bwt_line
  
#  fasta = path.basename(fasta)
#  fasta = fasta.replace('.fasta','')
#  fasta = fasta.replace('.fa','')
#  idxfile = "./idx/"+fasta+"."+name+"."+str(blocksize)
#  idxfile = fasta+"."+name+"."+str(blocksize)
#  np.save(idxfile+".IDX.npy",chr_index)
#  np.save(idxfile+".BWT.npy",chr_bwts)  
    #if n > 10: sys.exit(1)
  
