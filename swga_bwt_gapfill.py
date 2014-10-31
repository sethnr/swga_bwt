#!/usr/bin/python

import sys
sys.path.remove('/usr/lib/python2.7/dist-packages')
sys.path.append('/usr/lib/python2.7/dist-packages')


from Bio import SeqIO
from Bio.Seq import Seq
from math import floor,ceil,factorial
import random as rand
import numpy as np
import h5py as h5
import os.path as path
from re import match
import profile
import itertools
import pickle as pkl

import swga_bwt as s
import argparse

#############
# do some stuff
############

#backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
#           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']

parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
#parser.add_argument('-t|--target', action="store", dest='idxfile', type=str, help='target genome index', nargs='+')
parser.add_argument('-p|--patterns', action="store", dest='patternfile', type=str, help='list of patterns', nargs='+')
parser.add_argument('-i|--tmatchidx', action="store", dest='tmatchidx', type=str, help='match positions index', nargs='+')
parser.add_argument('-o|--out', action="store", dest='out', type=str, help='outfile', nargs='?')
parser.add_argument('-A|--array', action="store_true", dest='farm', help='run on farm? i.e. use LSB_JOBINDEX')

parser.add_argument('-c|--combs', action="store", dest='nCombs', default=1000, type=int, help='number of combinations to assess', nargs='?')
parser.add_argument('-n|--bestN', action="store", dest='bestN', default=100, type=int, help='show best N combinations', nargs='?')
parser.add_argument('-k|--plexSize', action="store", dest='plexSize', default=20, type=int, help='size of combinations to select', nargs='?')
parser.add_argument('-r|--minRatio', action="store", dest='minRatio', default=0, type=float, help='only include primers with higher than /r/ ratio in the plex', nargs='?')
parser.add_argument('-M|--noMaxPlex', action="store_false", dest='getMaxPlex', help='calculate max plex? (probably not worth doing on huge datasets')

args = parser.parse_args()

#pi = os.getenv("LSB_JOBINDEX")

farm = args.farm
#idxfile = args.idxfile[0]
#idxfile = idxfile.replace(".IDX.hdf5","")
patternfile = args.patternfile[0]
tmatchidx = args.tmatchidx[0]
getMaxPlex = args.getMaxPlex
nCombs = args.nCombs
bestN = args.bestN


k = args.plexSize
minRatio = args.minRatio # 10


#print >>sys.stderr, args.farm

if args.farm is True:
  pi = os.getenv("LSB_JOBINDEX")
  if pi == None:
    sys.exit(100,"pi = none - did you mean to submit an array job?")
  patternfile = patternfile+"."+pi

out = args.out
if out == None:
  out = sys.stdout
else:
  out = open(outfile,'w')


#	if len(sys.argv)==1: 
#	  print "\n\tpython swga_bwt_gapfill.py index_file pattern_file no_patterns_to_check:index\n\n"
#	  sys.exit(1)
#	
#	idxfile = sys.argv[1]
#	patternfile = sys.argv[2]
#	tmatchidx = sys.argv[3]
#	out = sys.stdout
#	
#	pblocksize=-1
#	if len(sys.argv) > 5:
#	  chr = sys.argv[4]
#	  pblocksize,pi = sys.argv[5].split(':')
#	  pblocksize = int(pblocksize)
#	  pi = int(pi)
#	elif len(sys.argv) > 4:
#	  pblocksize,pi = sys.argv[4].split(':')
#	  pblocksize = int(pblocksize)
#	  pi = int(pi)
#	#  out = open(outfile,'w')


#	pfile = open(patternfile,'r')
#	patterns = []
#	ratios = []
#	for line in pfile:
#	  if match('#',line):
#	    next
#	  else:
#	    F = line.split()
#	    if len(F) > 1:
#	      patterns += [F[0]]
#	      ratios += [F[-1]]
#	    else:
#	      patterns += [line.rstrip('\n').lower()]
#	      ratios += [1]

# limit patterns to pi'th set of pblocksize
#	if pblocksize != -1:
#	  pstart = pblocksize * pi
#	  pend = pstart + pblocksize
#	  if pstart >= len(patterns): sys.exit(1)
#	  elif pend >= len(patterns): pend = len(patterns)-1 
#	  patterns = patterns[pstart:pend]
#	
#	print >>sys.stderr, len(patterns)



      
#    print patterns
#ratios = np.array(ratios,dtype='float')
#patterns = ['nnnnn'] + patterns
#guess chrname and blocksize from indexname
# fasta, blocksize = path.basename(idxfile).split('.')

# blocksize = int(blocksize)
#chr_index = np.load(idxfile+".IDX.npy")
#chr_bwts = np.load(idxfile+".BWT.npy")

# index = h5.File(idxfile+".IDX.hdf5", "r")

#	if len(sys.argv) > 4:
#	  chrs = [chr]
#	else:
#	  chrs = index.keys()
#	
#	print >>sys.stdout, chrs


matchcounts = None
blockPosns = None

#	if path.exists(tmatchidx):
#	  print >> sys.stdout, "found file ",tmatchidx
#	  tmatchfile = open(tmatchidx,"r")
#	  (matchcounts, blockPosns) = pkl.load(tmatchfile)
#	  tmatchfile.close()
#	#  idxPatterns = sorted(matchcounts.keys(), key=lambda a: a[1])
#	  idxPatterns = set([key[1] for key in matchcounts.keys()])
#	  newPatterns = [p for p in patterns if p not in idxPatterns]
#	  print >>sys.stdout, idxPatterns,"\n",newPatterns, len(newPatterns)
#	  
#	  if(len(newPatterns) > 0):
#	    (newMatchCounts, newBlockPosns) = s.getMatchCounts(index,patterns,blocksize,chrs)
#	    matchcounts.update(newMatchCounts)
#	    cachecounts = open(tmatchidx, 'w')
#	    blockPosns.update(newBlockPosns)
#	    pkl.dump((matchcounts, blockPosns),cachecounts)
#	    cachecounts.close()
#	else:
#	  print >>sys.stdout, tmatchidx, "not found"
#	  (matchcounts, blockPosns) = s.getMatchCounts(index,patterns,blocksize,chrs)
#	  print >>sys.stdout, matchcounts
#	  cachecounts = open("match_counts.pkl","w")
#	  pkl.dump((matchcounts,blockPosns),cachecounts)
#	  cachecounts.close()

#	#initialise index if it doesn't exist
#	if not path.exists(tmatchidx):
#	  matchindex = h5.File(tmatchfile, "w")
#	#compression = None
#	  compression = 'gzip'
#	  matchcountsNP, blockPosnsNP, pArray = s.initMatchCountsHDF5(index,patterns,chrs)
#	  matchindex.close()
#	
#	
#	#update index with new patterns
#	tmatchfile = open(tmatchidx,"r")
#	matchindex = h5.File(tmatchfile, "a")
#	addMatchCountsHDF5(index, newPatterns, matchindex)
#	matchindex.close()
#	

#	  idxPatterns = set([key[1] for key in matchcounts.keys()])
#	  newPatterns = [p for p in patterns if p not in idxPatterns]
#	  print >>sys.stdout, idxPatterns,"\n",newPatterns, len(newPatterns)
#	  
#	  if(len(newPatterns) > 0):
#	    (newMatchCounts, newBlockPosns) = s.getMatchCounts(index,patterns,blocksize,chrs)
#	    matchcounts.update(newMatchCounts)
#	    cachecounts = open(tmatchidx, 'w')
#	    blockPosns.update(newBlockPosns)
#	    pkl.dump((matchcounts, blockPosns),cachecounts)
#	    cachecounts.close()
#	
#	  (matchcounts, blockPosns) = s.getMatchCounts(index,patterns,blocksize,chrs)
#	  print >>sys.stdout, matchcounts
#	  cachecounts = open("match_counts.pkl","w")
#	  pkl.dump((matchcounts,blockPosns),cachecounts)
#	  cachecounts.close()


#matchcounts, blockPosns = getMatchCounts(index,chrs)
#allBlocks = len(matchcounts.keys())
#print matchcounts
#matchcounts = np.recarray(matchcounts, dtype=[('pattern',"S15"),("count","i6")])
#print matchcounts.shape
#allBlocks = matchcounts.shape[0]



plex=[]
matchcountsHDF5 = h5.File(tmatchidx, "r")

allBlocks = matchcountsHDF5['matchcounts'].shape[0]
matchcounts = matchcountsHDF5['matchcounts']
blockPosns = matchcountsHDF5['blockposns']
patterns = matchcountsHDF5['patterns']
ratios = matchcountsHDF5['ratios']
emptyBlocks = np.array([True] * allBlocks)

print >>sys.stderr, getMaxPlex
if getMaxPlex == True:
#for all patterns, add pattern index to array if it improves on previous set
  print >>sys.stderr, "assessing ",len(patterns)," patterns for max plex"
  for i in range(0,len(patterns)):
#  print matchcounts[:,i].transpose()
    stillEmpty = np.sum(matchcounts[:,plex + [i]],axis=1)==0
    if i % 1000 == 0: print ".",
    if i % 20000 ==0: print ""
    if np.sum(stillEmpty) < np.sum(emptyBlocks):
      emptyBlocks = stillEmpty
      plex = plex + [i]
    if len(emptyBlocks)==0:
      break
else:
  plex = range(0,len(patterns))

filled = np.array(matchcounts > 0,dtype='bool')
sdcols = np.std(matchcounts,axis=0)
plex = np.array(plex,dtype='I5')
print >>sys.stderr, len(plex), filled.shape, ratios.shape, sdcols.shape

#stop=0

#combs = itertools.combinations(plex,k)
ratios = np.array(ratios)
patterns = np.array(patterns)
passRatio = np.greater_equal(ratios, [minRatio])

#print >>sys.stderr, sum(ratios[plex] > minRatio)
print >>sys.stderr, passRatio.shape
#plex = plex[[ratios[plex] >= minRatio]]
plex = plex[passRatio]

n = len(plex)
print n, k, n-k
#print n, factorial(n)
#print k, factorial(k)
#print n-k, factorial(n-k)

noCombs = factorial(n) / (factorial(k) * factorial(n-k))
print >> sys.stderr, noCombs, " total combinations of ",k," in ",n," elements"
#combs = tuple(combs)

#for comb in combs:
#len(plex)

#make normal probability distribution from first in pattern to last
mu, sigma = 0.0, 0.1 # mean and sd
probs = np.sort(abs(np.random.normal(loc=mu, scale=sigma, size=len(plex))))[::-1]
probs = probs/sum(probs)

print >> sys.stderr, "asessing ", nCombs, "random samples"

ranks = {}
i = 0
while i <= nCombs:
#  comb = rand.sample(plex,k)
  comb = np.sort(np.random.choice(plex, k, replace=False, p=probs))
  i += 1
#  if i > 5000: break
  if i % 50 == 0: print >> sys.stderr, '.',
  if i % 1000 == 0: print >> sys.stderr, i," of ",noCombs
#  filledComb = filled[:,comb]
#  filledRow = np.apply_along_axis(any,arr=filledComb,axis=1)
  filledRow = np.sum(matchcounts[:,comb],axis=1)>0
#  print sum(filledRow), sum(filledRowSum)
  sdcomb = sdcols[np.array(comb)]

  print >>sys.stderr, np.array(comb)
  
  meanratio = np.mean(ratios[np.array(comb)])
  
  if 0 == 1:
    print "filled[:,comb]"
    profile.run("filled[:,comb]")
    print "np.apply_along_axis(any,arr=filledComb,axis=1)".upper()
    profile.run("np.apply_along_axis(any,arr=filled[:,comb],axis=1)")
    print "ROWSUMMING"
    profile.run("np.sum(matchcounts[:,comb],axis=1)>0")

  unfilled = len(filledRow) - sum(filledRow)
  #divide by ten to get rough number of unfilled blocks
  ratioBlock = floor(1/meanratio*1000)/1000
  key = (floor(unfilled/10),ratioBlock,np.mean(sdcomb), i)
  ranks[key] = (comb, unfilled, meanratio)

  
ranks.keys()

bestorder = sorted(ranks.keys(), key=lambda a: a[0:2])

i = 0;
for key in bestorder:
  i +=1
  print key
  plex, something, somethingelse = ranks[key]
  
  print "ranks ", ranks[key]
  print "key", key
  print "plex ", plex
  print "\t".join(list(key)) + "\t",
  print join(patterns[plex])
  if i == bestN: break
sys.exit(0)
