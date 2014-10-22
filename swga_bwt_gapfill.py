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

import swga_bwt


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

if len(sys.argv) > 3:
  chrs = [chr]
else:
  chrs = index.keys()


if path.exists(tmatchidx):

  tmatchfile = open(tmatchidx,"w")
  matchcounts, blockPosns, pArray = pkl.load(tmatchidx)
  tmatchfile.close()
  primIndices = list(matchcounts.keys(), key=lambda a: a[1])

  #bestorder = sorted(ranks.keys(), key=lambda a: a[0:2])
  ixdPatterns = allPatterns[primIndices]
  newPatterns = [p for p in patterns if p not in idxPatterns]
  if(len(newPatterns) >= 0):
    newMatchCounts, newBlockPosns = getMatchCounts(index,patterns,blocksize,chrs)
    matchcounts.update(newMatchCounts)
    cachecounts = open(tmatchfile, 'w')
    pkl.dump((matchcounts, block_posns),cachecounts)
    cachecounts.close()
else:
  matchcounts, blockPosns = getMatchCounts(index,patterns,blocksize,chrs)
  cachecounts = open("match_counts.pkl","w")
  pkl.dump((matchcounts,blockPosns),cachecounts)
  cachecounts.close()



#matchcounts, blockPosns = getMatchCounts(index,chrs)
allBlocks = matchcounts.shape[0]
plex=[]
emptyBlocks = np.array([True] * allBlocks)

for i in range(0,len(patterns)):
#  print matchcounts[:,i].transpose()
  stillEmpty = np.sum(matchcounts[:,plex + [i]],axis=1)==0

#  print stillEmpty
  if np.sum(stillEmpty) < np.sum(emptyBlocks):
#    print >> sys.stderr, i, np.sum(emptyBlocks), np.sum(stillEmpty), len(plex),
#    print >> sys.stderr, plex
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

#  print ratios[:10]
#  print ratios[np.array(comb)]
#  print np.mean(ratios[np.array(comb)])
  meanratio = np.mean(ratios[np.array(comb)])
  
  if 0 == 1:
    print "filled[:,comb]"
    profile.run("filled[:,comb]")
    print "np.apply_along_axis(any,arr=filledComb,axis=1)".upper()
    profile.run("np.apply_along_axis(any,arr=filled[:,comb],axis=1)")
    print "ROWSUMMING"
    profile.run("np.sum(matchcounts[:,comb],axis=1)>0")
#  print "sdcols[np.array(comb)]"
#  profile.run("sdcols[np.array(comb)]")
#  print "np.mean(ratios[np.array(comb)])"
#  profile.run("np.mean(ratios[np.array(comb)])")
#  print sdcols[:20]
#  print comb
#  print sdcomb
#  print np.mean(sdcomb)
  unfilled = len(filledRow) - sum(filledRow)
  #divide by ten to get rough number of unfilled blocks
  ratioBlock = floor(1/meanratio*1000)/1000
  key = (floor(unfilled/10),ratioBlock,np.mean(sdcomb), i)
  ranks[key] = (comb, unfilled, meanratio)
#  print ranks[key], key
  
ranks.keys()
#print ranks.values()
bestorder = sorted(ranks.keys(), key=lambda a: a[0:2])

i = 0;
for key in bestorder:
  i +=1
  print key, ranks[key]
  if i == 100: break
sys.exit(1)
