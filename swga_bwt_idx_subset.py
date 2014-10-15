#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
import numpy as np
import h5py as h5
import os.path as path
import gzip
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
percent = float(sys.argv[3])

print fasta

if fasta[-2:] == "gz": 
  seqfile = gzip.open(fasta)
else:
  seqfile = open(fasta)

seq = SeqIO.parse(seqfile,'fasta')

idxfile = path.basename(fasta)
idxfile = idxfile.replace('.fasta','')
idxfile = idxfile.replace('.fa','')
idxfile = idxfile.replace('.gz','')
idxfile = "./idx/"+idxfile+"."+str(percent)+"."+str(blocksize)

index = h5.File(idxfile+".IDX.hdf5", "a")
compression = 'gzip'

#ci = 0
genomeSize=0
genomeBlocks=0

for chr in seq:
  seqlen = len(chr.seq)
  print chr.name, seqlen
  genomeSize += seqlen
  lastBlock = (seqlen/blocksize)
  genomeBlocks += lastBlock
totalBlocks = int((genomeBlocks*percent))

print "blocks",totalBlocks, blocksize+1
print genomeSize*percent, int(floor(genomeBlocks*percent))

chr_index = index.create_dataset("subset/idx",
                                 (totalBlocks+1,blocksize+1,len(allBases)), 
                                 dtype='i',
                                 compression=compression)
chr_bwts = index.create_dataset("subset/bwt", 
                                (totalBlocks+1,blocksize+1), 
                                dtype='S1',
                                compression=compression)

seqfile.seek(0) #reset file pointer
seq = SeqIO.parse(seqfile,'fasta')
bi=0  #blockindex
for chr in seq:
  ci=0
  name = chr.id
  seqlen = len(chr.seq)
#  print name, seqlen
  lastBlock = (seqlen/blocksize)
    #lastBlock #rounds down as both are integers

  for n in range(0,lastBlock):
    bi +=1

    if bi % (1/percent) == 0:
      end = (n+1)*blocksize
      ci +=1
      if end > seqlen: #if not full block, pad with Ns
        pad = 'n'*(end-seqlen)
        seqblock = chr.seq[n*blocksize:seqlen]
        seqblock = seqblock + Seq(pad)
      else:  
        seqblock = chr.seq[n*blocksize:end]
      
#    print n, n*blocksize, end, len(seqblock)
      bwt_line = bwt(str(seqblock.lower()))     
      baseRanks = rankAllBwtNP(bwt_line)
      i = int(bi/(1/percent))
      chr_index[i] = baseRanks
      chr_bwts[i] = bwt_line
 # print name, ci, "blocks indexed"
print int(bi/(1/percent)), "blocks indexed"

#  fasta = path.basename(fasta)
#  fasta = fasta.replace('.fasta','')
#  fasta = fasta.replace('.fa','')
#  idxfile = "./idx/"+fasta+"."+name+"."+str(blocksize)
#  idxfile = fasta+"."+name+"."+str(blocksize)
#  np.save(idxfile+".IDX.npy",chr_index)
#  np.save(idxfile+".BWT.npy",chr_bwts)  
    #if n > 10: sys.exit(1)
  
