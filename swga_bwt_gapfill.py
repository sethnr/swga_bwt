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


#############
# do some stuff
############

#backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
#           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']


if len(sys.argv)==1: 
  print "\n\tpython swga_bwt_gapfill.py index_file pattern_file no_patterns_to_check:index\n\n"
  sys.exit(1)

idxfile = sys.argv[1]
patternfile = sys.argv[2]
tmatchidx = sys.argv[3]
out = sys.stdout

pblocksize=-1
if len(sys.argv) > 5:
  chr = sys.argv[4]
  pblocksize,pi = sys.argv[5].split(':')
  pblocksize = int(pblocksize)
  pi = int(pi)
elif len(sys.argv) > 4:
  pblocksize,pi = sys.argv[4].split(':')
  pblocksize = int(pblocksize)
  pi = int(pi)
#  out = open(outfile,'w')


pfile = open(patternfile,'r')
patterns = []
ratios = []
for line in pfile:
  if match('#',line):
    next
  else:
    F = line.split()
    if len(F) > 1:
      patterns += [F[0]]
      ratios += [F[-1]]
    else:
      patterns += [line.rstrip('\n').lower()]
      ratios += [1]

# limit patterns to pi'th set of pblocksize
if pblocksize != -1:
  pstart = pblocksize * pi
  pend = pstart + pblocksize
  if pstart >= len(patterns): sys.exit(1)
  elif pend >= len(patterns): pend = len(patterns)-1 
  patterns = patterns[pstart:pend]

print >>sys.stderr, len(patterns)



      
#    print patterns
ratios = np.array(ratios,dtype='float')
#patterns = ['nnnnn'] + patterns
#guess chrname and blocksize from indexname
fasta, blocksize = path.basename(idxfile).split('.')

blocksize = int(blocksize)
#chr_index = np.load(idxfile+".IDX.npy")
#chr_bwts = np.load(idxfile+".BWT.npy")

index = h5.File(idxfile+".IDX.hdf5", "r")

if len(sys.argv) > 4:
  chrs = [chr]
else:
  chrs = index.keys()

print >>sys.stdout, chrs


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

#initialise index if it doesn't exist
if not path.exists(tmatchidx):
  matchindex = h5.File(tmatchfile, "w")
#compression = None
  compression = 'gzip'
  matchcountsNP, blockPosnsNP, pArray = s.initMatchCountsHDF5(index,patterns,chrs)
  matchindex.close()


#update index with new patterns
tmatchfile = open(tmatchidx,"r")
matchindex = h5.File(tmatchfile, "a")
addMatchCountsHDF5(index, newPatterns, matchindex)
matchindex.close()


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

allBlocks = matchcountsHDF5['matchcounts'].shape[0]
matchcounts = matchcountsHDF5['matchcounts']
blockPosns = matchcountsHDF5['block_posns']
allPatterns = matchcountsHDF5['patterns']

emptyBlocks = np.array([True] * allBlocks)


#for all patterns, add pattern index to array if it improves on previous set
for i in range(0,len(patterns)):
#  print matchcounts[:,i].transpose()
  stillEmpty = np.sum(matchcounts[:,plex + [i]],axis=1)==0
  if np.sum(stillEmpty) < np.sum(emptyBlocks):
    emptyBlocks = stillEmpty
    plex = plex + [i]
  if len(emptyBlocks)==0:
    break


filled = np.array(matchcounts > 0,dtype='bool')
sdcols = np.std(matchcounts,axis=0)
plex = np.array(plex,dtype='I5')
print len(plex),filled.shape,sdcols.shape

stop=0
ranks = {}
i = 0
k = 30
minRatio = 0 # 10
combs = itertools.combinations(plex,k)
print sum(ratios[plex] > minRatio)
plex = plex[[ratios[plex] >= minRatio]]


n = len(plex)
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
#print probs
noSamples = 10000
print >> sys.stderr, "asessing ", noSamples, "random samples"

while i <= noSamples:
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
  print key, ranks[key]
  if i == 100: break
sys.exit(1)
