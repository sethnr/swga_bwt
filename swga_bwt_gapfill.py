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

allBases = ['a','c','g','t','n','$']
allBases.sort()
# nb: ^^ allBases MUST be sorted - next in array must be next in bwt matrix


def firstColNP(tots):
  # get ranges of each base in first col (contiguous as BWM is sorted)
  first = {}
  totc = 0
  for n in range(0,len(allBases)):
    c = allBases[n]
    count = tots[n]
    first[c] = (totc, totc + count)
    totc += count
  return first

def countMatchesNP(baseRanks, first, p):
  if p[-1] not in first:
    return 0 # character doesn't occur in T
  l, r = first[p[-1]]
  i = len(p)-2
  while i >= 0 and r > l:
    c = p[i]
    ci = allBases.index(c)
    l = first[c][0] + baseRanks[l-1,ci]
    r = first[c][0] + baseRanks[r-1,ci]
    i -= 1
  return r - l # return size of final range


#############
# do some stuff
############

backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']

idxfile = sys.argv[1]
patternfile = sys.argv[2]

out = sys.stdout

if len(sys.argv) > 3:
  chr = sys.argv[3]
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
  chrs = sys.argv[3:]
else:
  chrs = index.keys()



#for b in backgrounds:
#  background = h5file(b,"r")
#  b_index = index["subset/idx"]
#  b_bwts = index["subset/bwt"]
#  for p in patterns:
#    

def getMatchCounts(index, chrs):
  allBlocks = 0
  for chr in chrs:
    chr_bwts = index[chr+"/bwt"]
    allBlocks += chr_bwts.shape[0]

  matchcounts = np.zeros(shape=(allBlocks,len(patterns)),dtype='int')
  blockPosns = np.recarray(shape=(allBlocks,3),dtype=[('chr', 'S12'), ('start', int),('end',int)])

#def getMatchCounts(index, chrs):
  #sys.exit(1)
  bi=-1
  for chr in chrs:
    chr_index = index[chr+"/idx"]
    chr_bwts = index[chr+"/bwt"]

    print >> sys.stderr, fasta, blocksize, chr
    blocks = chr_index.shape[0]

    for n in range(0,blocks):
      baseRanks = chr_index[n]
      bwt_line = chr_bwts[n]
      bi +=1
    #get mapping of each base to range of posns in first column 
    #(all contiguous as col is sorted)
      firstColMap = firstColNP(baseRanks[-1])
    
    #    print >> out, chr, n*blocksize, (n+1)*blocksize, blocksize,
      pi = -1
      blockPosns[bi,0:2] = (chr,n*blocksize,(n+1)*blocksize)
      for p in patterns:
        pi +=1
        matchcounts[bi,pi] = countMatchesNP(baseRanks,firstColMap,p)
#      print >> out, countMatchesNP(baseRanks,firstColMap,p),
  return matchcounts, blockPosns
#    print >> out, ""


if path.exists("match_counts.pkl"):
  countfile = open("match_counts.pkl","r")
  (matchcounts, blockPosns) = pkl.load(countfile)
  countfile.close()
else:
  matchcounts, blockPosns = getMatchCounts(index,chrs)
  countfile = open("match_counts.pkl","w")
  pkl.dump((matchcounts,blockPosns),countfile)
  countfile.close()

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
k = 10
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
print >> sys.stderr, "asessing ", 5000, "random samples"

#for comb in combs:
#len(plex)
mu, sigma = 0.0, 0.1 # mean and sd
probs = np.sort(abs(np.random.normal(loc=mu, scale=sigma, size=len(plex))))[::-1]
probs = probs/sum(probs)
#print probs
noSamples = 5000

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
