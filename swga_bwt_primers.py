#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
from itertools import permutations, product
import numpy as np
import h5py as h5
import os.path as path
import gzip
from re import findall
#import profile

bases = ['a','c','t','g']
#print itertools.permutations(bases,8)

#primerLen = 8
#primers = list(map("".join, product(bases, repeat=primerLen)))
#print len(primers)

tm_limit = 25
lower_len=12
upper_len=12

target = './idx/Pf3D7_v3.0.01.3000'
backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.0.01.3000',
           './idx/Homo_sapiens.GRCh38.dna_sm.primary_assembly.0.001.3000']

#idxfile = sys.argv[1]
#patternfile = sys.argv[2]
filterLength=False

seqfile = sys.argv[1]


def getTM(primer):
    nA = len(findall('a',primer))
    nC = len(findall('c',primer))
    nT = len(findall('t',primer))
    nG = len(findall('g',primer))
    tm = nA*2+nT*2+nC*4+nG*4
    #print tm
    return tm


def getIndexCounts(target, backgrounds, patterns):
  backcounts = {}

  for b in [target] + background:
    b = +".IDX.hdf5"
    print b
    bgindex = h5.File(b,"r")
    b_index = bgindex["subset/idx"]
    b_bwts = bgindex["subset/bwt"]
  #firstColMap = firstColNP(baseRanks[-1])
    for p in patterns:
    #print >> out, chr, n*blocksize, (n+1)*blocksize, blocksize,
      noBlocks = b_bwts.shape[0]
      noBlocks = 10
      matches = 0
      print p,
      for n in range(0,noBlocks):
        baseRanks = b_index[n]
        firstColMap = firstColNP(baseRanks[-1])
        bwt_line = b_bwts[n]
        match = countMatchesNP(baseRanks,firstColMap,p)
        print match,
        matches += match
        backcounts[p,b] = float(matches) / float(noBlocks)
      print ''


seq = SeqIO.parse(seqfile,'fasta')

primers = set()
for chr in seq:
    seqstr = str(chr.seq).lower()
    seqlen = len(seqstr)
    print chr.name, seqlen,
    for i in range(0,seqlen):
        for primerLen in reversed(range(lower_len,upper_len+1)):
            primer = seqstr[i:(i+primerLen)]
#        print i, i+primerLen, len(primer), primer
            if len(primer) == primerLen:
                tm = getTM(primer)
                if tm <= tm_limit:
                    primers.add(primer)
    print len(primers)


if not filterLength:
    print "no length filtration"
else:
    primersF = set(primers)
    for primer1 in primers:
        found = False
        for primer2 in primersF:
            if primer2 != primer1 and primer2 in primer1:
                found = True
                primersF.remove(primer2)
#                print primer1, "contains", primer2
            if found: break
    primers = primersF
print len(primers)

out = open("candidates.txt","w")
print >> out, "\n".join(primers)



#backcounts = getIndexCounts(target, background, patterns)

sys.exit(1)

