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
import swga_bwt as s
import argparse

allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix

########
# get some vars
########
parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-f','--fasta', action="store", dest='fasta', type=str, help='target genome to index (fasta)', nargs='?')
parser.add_argument('-o','--out', action="store", dest='out', type=str, help='outfile', nargs='?')
parser.add_argument('-I','--idxdir', action="store", dest='idxdir', type=str, help='directory for indices', nargs='?', default="./idx/")
parser.add_argument('-b','--blocksize', action="store", dest='blocksize', type=int, default=3000, help='block size (kb)', nargs='?')
parser.add_argument('-p','--background_percent', action="store", dest='bpc', type=float, help='percent of genome to index', default=100, nargs='?')

args = parser.parse_args()


fasta = args.fasta
blocksize=args.blocksize
percent = float(args.bpc)/100
idxdir = args.idxdir

#############
# do some stuff
############

#blocksize = 5000

# fasta = sys.argv[1]
# blocksize = int(sys.argv[2])
# percent = float(sys.argv[3])

print 1/percent

if fasta[-2:] == "gz": 
  seqfile = gzip.open(fasta)
else:
  seqfile = open(fasta)

seq = SeqIO.parse(seqfile,'fasta')

idxfile = path.basename(fasta)
idxfile = idxfile.replace('.fasta','')
idxfile = idxfile.replace('.fa','')
idxfile = idxfile.replace('.gz','')
idxfile = idxdir+idxfile+"."+str(percent)+"."+str(blocksize)

print >>sys.stderr, idxfile
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
      bwt_line = s.bwt(str(seqblock.lower()))     

#      print >>sys.stderr,seqblock.lower(),"\n",bwt_line
        
      baseRanks = s.rankAllBwtNP(bwt_line)
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
  
