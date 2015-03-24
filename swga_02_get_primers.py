#!/usr/bin/python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil
#from itertools import permutations, product #, iter
import itertools as it
import numpy as np
import h5py as h5
import os.path as path
import gzip
# from re import findall
#import profile
import argparse
from swga_bwt import * 
from myParallel import *

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return it.izip_longest(fillvalue=fillvalue, *args)

parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
# parser.add_argument('-t','--target', action="store", dest='idxfile', type=str, help='target genome index', nargs='?')
parser.add_argument('-b','--background', action="store", dest='backgrounds', type=str, help='background genomes', nargs='+')
parser.add_argument('-f','-t','--fasta', action="store", dest='fasta', type=str, help='target genome (fasta)', nargs='?')
# parser.add_argument('-p|--patterns', action="store", dest='patternfile', type=str, help='list of patterns', nargs='?')
# parser.add_argument('-A|--array', action="store_true", dest='farm', help='run on farm? i.e. use LSB_JOBINDEX')

parser.add_argument('-n|--minLength', action="store", dest='lower_len', default=8, type=int, help='lower size limit', nargs='?')
parser.add_argument('-x|--maxLength', action="store", dest='upper_len', default=12, type=int, help='lower size limit', nargs='?')
parser.add_argument('-c|--maxTm', action="store", dest='tm_limit', default=30, type=int, help='only include primers withlower than \'c\' melting temperature in the plex (c)', nargs='?')

#parser.add_argument('-r|--minRatio', action="store", dest='minRatio', default=10, type=float, help='only include primers with higher than /r/ ratio in the plex', nargs='?')

parser.add_argument('--lengthFilter', action="store", dest='filterLength', default=False, type=bool, help='filter primers for those contained within another (prioritise longer)', nargs='?')

parser.add_argument('-o','--out', action="store", dest='outfile', type=str, default="candidates.txt", help='target genome fasta', nargs='?')
parser.add_argument('-B','--blocksize', action="store", dest='blocksize', default=100, type=int, help='size of blocks for threading', nargs='?')
parser.add_argument('-T','--threads', action="store", dest='threads', default=10, type=int, help='no of threads [4]', nargs='?')

args = parser.parse_args()

#farm = args.farm
outfile = args.outfile


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
#target = args.idxfile
threads = args.threads
backgrounds = args.backgrounds

filterLength=args.filterLength

seqfile = args.fasta
blocksize = args.blocksize
blockToRPK = float(1000) / blocksize

seq = SeqIO.parse(seqfile,'fasta')

primers = set()
tcount = {}
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
                        tcount[primer] += 1
                    else:
                        primers.add(primer)
                        tcount[primer] = 1
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


out = open(outfile,"w")

#print >>out, "#primer","total","rpk"
#for primer in primers:
#    print >>out, primer, tcount[primer], (float(tcount[primer]) / genomeLen) * 1000


#def printIndexCountsBlock(primerBlock, outfile):
#    out = open(outfile,'w')
#    for primer in block:
#        backcounts = getIndexCounts(target, background, patterns)
#        print >>out, primer, tcount[primer], (float(tcount[primer]) / genomeLen) * 1000


#use N threads to get counts in backgrounds 
backcounts = {}


def _getMatchesBlock(newPatterns, *args):
    global backcounts
    print >>sys.stderr, len(newPatterns)
    if(len(newPatterns) > 0):
        newCounts = getIndexCounts(backgrounds,newPatterns, blockToRPK)
        backcounts.update(newCounts)
    
if threads > 1:
    primerBlocks = grouper(primers, blocksize)

<<<<<<< HEAD
    print >>sys.stderr, blocksize, primerBlocks
    primerBlocksL = list(primerBlocks)
=======
  #  print >>sys.stderr, blocksize, primerBlocks
  #  primerBlocksL = list(primerBlocks)
  #  for i in range(1,5):
  #      print >>sys.stderr, primerBlocksL[i]
>>>>>>> 76bf2afd67e968536cb7c15e9602121f98b8cb4d
    parallel(_getMatchesBlock, threads, primerBlocks)
else:
    print >>out, "#primer","total","rpk"
    for primer in primers:
        print >>out, primer, tcount[primer], (float(tcount[primer]) / genomeLen) * 1000



        
for p in primers:
#  print p, target
#  tcount = backcounts[p,target+".IDX.hdf5"]
  tcount = tcount[p]
  ratioTotal = 0
  print >>out, p, tcount, 
  #for b in backgrounds.values():
  for b in backgrounds:
    bcount = backcounts[(p,b+".IDX.hdf5")]
    if bcount > 0: ratio = tcount/bcount
    else: ratio=float('inf')
    print >>out, bcount, ratio,
    ratioTotal += ratio
  print >>out, ratioTotal / len(backgrounds)
#  print >>out, ''

        
#print >> out, "#primer", "t_count",
#for b in backgrounds.keys():
#  print >>out, str(b)+"_count", str(b)+"_ratio",
#print >> out, "mean_ratio"

        

#print >> out, "\n".join(primers)

#backcounts = getIndexCounts(target, background, patterns)

sys.exit(0)

