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
import argparse


parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
parser.add_argument('-t|--target', action="store", dest='idxfile', type=str, help='target genome index', nargs='?')
parser.add_argument('-b|--backgrounds', action="store", dest='backgrounds', type=str, help='list of patterns', nargs='+')
parser.add_argument('-f|--fasta', action="store", dest='fasta', type=str, help='target genome fasta', nargs='?')
parser.add_argument('-p|--patterns', action="store", dest='patternfile', type=str, help='list of patterns', nargs='?')
parser.add_argument('-A|--array', action="store_true", dest='farm', help='run on farm? i.e. use LSB_JOBINDEX')

parser.add_argument('-n|--minLength', action="store", dest='lower_len', default=8, type=int, help='lower size limit', nargs='?')
parser.add_argument('-x|--maxLength', action="store", dest='upper_len', default=12, type=int, help='lower size limit', nargs='?')
parser.add_argument('-c|--maxTm', action="store", dest='tm_limit', default=30, type=int, help='only include primers withlower than \'c\' melting temperature in the plex (c)', nargs='?')

parser.add_argument('-r|--minRatio', action="store", dest='minRatio', default=10, type=float, help='only include primers with higher than /r/ ratio in the plex', nargs='?')

parser.add_argument('--lengthFilter', action="store", dest='filterLength', default=False, type=boolean, help='filter primers for those contained within another (prioritise longer)', nargs='?')

args = parser.parse_args()

farm = args.farm



bases = ['a','c','t','g']

#tm_limit = 30
#lower_len=8
#upper_len=12

#target = './idx/Pf3D7_v3.0.01.3000'
#backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.0.01.3000',
#           './idx/Homo_sapiens.GRCh38.dna_sm.primary_assembly.0.001.3000']
tm_limit = args.tm_limit
lower_len = args.lower_len
upper_len = args.upper_len
target = args.idxfile
backgrounds = args.backgrounds

filterLength=args.filterLength

seqfile = args.fasta


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
pcount = {}
genomeLen = 0
for chr in seq:
    seqstr = str(chr.seq).lower()
    seqlen = len(seqstr)
    print chr.name, seqlen,
    genomeLen += seqlen
    for i in range(0,seqlen):
        for primerLen in reversed(range(lower_len,upper_len+1)):
            primer = seqstr[i:(i+primerLen)]
#        print i, i+primerLen, len(primer), primer
            if len(primer) == primerLen:
                tm = getTM(primer)
                if tm <= tm_limit:
                    if primer in primers:
                        pcount[primer] += 1
                    else:
                        primers.add(primer)
                        pcount[primer] = 1
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

print >>out, "#primer","total","rpk"
for primer in primers:
    print >>out, primer, pcount[primer], (float(pcount[primer]) / genomeLen) * 1000

#print >> out, "\n".join(primers)

#backcounts = getIndexCounts(target, background, patterns)

sys.exit(1)

