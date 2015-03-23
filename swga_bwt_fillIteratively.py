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
import os
import swga_bwt as s
import argparse
# from sklearn.decomposition import PCA
# import matplotlib.pyplot as plt


#############
# do some stuff
############

#backgrounds = ['./idx/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.3000.IDX.hdf5',
#           'Homo_sapiens.GRCh38.dna_sm.primary_assembly.3000.IDX.hdf5']

parser = argparse.ArgumentParser(description='build index of genome matches in 3k blocks')
#parser.add_argument('-t|--target', action="store", dest='idxfile', type=str, help='target genome index', nargs='+')
parser.add_argument('-p|--patterns', action="store", dest='patternfile', type=str, help='list of patterns', nargs='+')
parser.add_argument('-i|--tmatchidx', action="store", dest='tmatchidx', type=str, help='match positions index', nargs='+')
parser.add_argument('-o|--out', action="store", dest='out', type=str, help='outfile', nargs='?',default="outfile")
parser.add_argument('-A|--array', action="store_true", dest='farm', help='run on farm? i.e. use LSB_JOBINDEX')

parser.add_argument('-c|--combs', action="store", dest='nCombs', default=1000, type=int, help='number of combinations to assess', nargs='?')
parser.add_argument('-n|--bestN', action="store", dest='bestN', default=100, type=int, help='show best N combinations', nargs='?')
parser.add_argument('-k|--plexSize', action="store", dest='plexSize', default=20, type=int, help='size of combinations to select', nargs='?')
parser.add_argument('-r|--minRatio', action="store", dest='minRatio', default=10, type=float, help='only include primers with higher than /r/ ratio in the plex', nargs='?')
parser.add_argument('-M|--noMaxPlex', action="store_false", dest='getMaxPlex', help='calculate max plex? (probably not worth doing on huge datasets')
parser.add_argument('-F|--forwardWeight', action="store", dest='forwardWeight', type=int, default = 4, help='weighting towards first primers (higher human/plasmo frequencies')
parser.add_argument('-C|--minCount', action="store", dest='minCount', default=100, type=int, help='only include primers with higher than /C/ counts in the target genome', nargs='?')
parser.add_argument('-D|--decSteps', action="store", dest='decs', default=5, type=int, help='decrease count and ratio by 1/D for d steps when chosing suboptimal probes', nargs='?')


args = parser.parse_args()

#pi = os.getenv("LSB_JOBINDEX")

farm = args.farm
#idxfile = args.idxfile[0]
#idxfile = idxfile.replace(".IDX.hdf5","")
#patternfile = args.patternfile[0]
tmatchidx = args.tmatchidx[0]
getMaxPlex = args.getMaxPlex
nCombs = args.nCombs
bestN = args.bestN
decs = args.decs
forwardWeight = args.forwardWeight
outfile = args.out

k = args.plexSize
minRatio = args.minRatio # 10
minCount = args.minCount # 10

if minRatio % 1 == 0:
  outfile = outfile+"_k"+str(k)+"_r"+str(int(minRatio))+"_c"+str(minCount)+"_d"+str(decs)
else:
  outfile = outfile+"_k"+str(k)+"_r"+str(minRatio)+"_c"+str(minCount)+"_d"+str(decs)


matchcounts = None
blockPosns = None

print >>sys.stderr, tmatchidx
#plex=[]
matchcountsHDF5 = h5.File(tmatchidx, "r")

allBlocks = matchcountsHDF5['matchcounts'].shape[0]
matchcounts = matchcountsHDF5['matchcounts']
blockPosns = matchcountsHDF5['blockposns']
blockchrs = matchcountsHDF5['blockchrs']
patterns = matchcountsHDF5['patterns']
ratios = np.array(matchcountsHDF5['ratios'])
emptyBlocks = np.array([True] * allBlocks)

filled = np.greater(matchcounts[:,:],[0])
sdcols = np.std(matchcounts,axis=0)

passFilters = []

print >>sys.stderr, filled.shape, ratios.shape, sdcols.shape

logfile = open(outfile+".log",'w')

# internal method to reset matrix - assumes filled/matchcounts/etc already there
# returns matrix like filled but with unfillable blocks set to false 
# and already filled blocks set to false
def _setMatrix(minRatio, minCount, plex):
  global passFilters

#  print >>sys.stderr, "resetting filled matrix"
  testFilled = np.copy(filled)

  passRatio = np.greater_equal(ratios, [minRatio])
  failRatio = np.less(ratios, [minRatio])
#  print "failRatio:",sum(failRatio)
  testFilled[:,failRatio] = False #remove those which are too low ratio

  patternCounts = np.array(matchcounts).sum(axis=0)
  passCount = np.greater_equal(patternCounts, [minCount])
  failCount = np.less(patternCounts, [minCount])
  testFilled[:,failCount] = False #remove those which are too low count
#  print "failCount:",sum(failCount)

  passFilters = np.vstack([passRatio,passCount]).all(axis=0)
#  print "passFilters:",sum(passFilters)

  indices = [p[-1] for p in plex]
  if (len(indices) > 0):
    #print >> sys.stderr, indices
    testFilled[:,np.array(indices)] = False
  #print >>sys.stderr, minRatio, minCount, plex, sum(sum(testFilled))
  return testFilled

def _getBestN(testFilled, filledPlex, k, reset=True):
#  global filledBlocksPotential
  plex = []
  pi = 0
  while pi < k:
    filledTot = testFilled.sum(axis=0)
    nextBestI = filledTot.argsort()[-1]
#  filledPlex
#  filled[:,nextBestI]
    # set newly filled to true in 'filled' array
    filledPlex[filled[:,nextBestI]] = True
    filledBlocksPotential =  _getTotalFills(passFilters,filled)

    filledBlocksPrimer = _getTotalFills([nextBestI],filled)
    #blank out array from blocks yet to fill, then calc how many are left
    testFilled[filledPlex,:] = False
    filledBlocksRemain =  _getTotalFills(passFilters,testFilled)
    #add new line to plex
    pLine = (patterns[nextBestI],
       ratios[nextBestI],
       patternCounts[nextBestI],
       filledBlocksPrimer,
       filledBlocksPotential,
       nextBestI)
    plex.append(pLine)
    pi += 1
    
    filledBlocksPlex = _getTotalFills([p[-1] for p in plex],filled)

    print pi, str(filledBlocksPlex)+"/"+str(filledBlocksRemain), str(filledBlocksPotential)+"/"+str(allBlocks), 
    print patternCounts[nextBestI], patterns[nextBestI]
#    print filledBlocksPlex, sum(filledPlex), filledPlex.shape
#    if sum(filledPlex) == filledBlocksPotential:
    if filledBlocksRemain == 0:
      print >>sys.stderr, "no more fillable gaps"
      if reset:
        filledPlex = np.zeros((allBlocks,),dtype=bool) 
        testFilled = _setMatrix(minRatio, minCount, plex)
#    print >>sys.stderr, sum(filledPlex), allBlocks
      else:
        return plex
  return plex

def _getTotalFills(plex, thisFilled):
#  plexi = np.array([p[-1] for p in plex])
#  print >>sys.stderr, plex
  plexi = np.array(plex)
  filled = sum(thisFilled[:,plexi].any(axis=1))
  return filled


#get 1D array of 'is block filled' (all false)
filledPlex = np.zeros((allBlocks,),dtype=bool) 

testFilled = _setMatrix(minRatio, minCount, [])
patternCounts = np.array(matchcounts).sum(axis=0)

#get optimal probes:
filledBlocksPotential = _getTotalFills(passFilters,filled)

plex=_getBestN(testFilled,filledPlex,k)

filledPlexCount = _getTotalFills([p[-1] for p in plex],filled)
plexRatios = ratios[np.array([p[-1] for p in plex])]
meanRatioOpt = plexRatios[np.isfinite(plexRatios)].mean()

print >>logfile, outfile, filledBlocksPotential, filledPlexCount, (allBlocks - filledPlexCount), meanRatioOpt, 


#get suboptimal probes

decrementRatio = minRatio/decs
decrementCount = float(minCount)/decs
subopt = []
n = 1
#reset filled matrix to remove already-chosen probes
testFilled = _setMatrix((minRatio - n*decrementRatio), 
                          (minCount - n*decrementCount), 
                          plex)

#print >>sys.stderr, sum(filledPlex), allBlocks
while _getTotalFills([p[-1] for p in (plex + subopt)],filled) < allBlocks:
  suboptRatio = (minRatio - n*decrementRatio)
  suboptCount = (minCount - n*decrementCount)
  print >>sys.stderr, "getting suboptimal probes r",suboptRatio," c",suboptCount
  print sum(sum(testFilled)), "->",
  testFilled[filledPlex,:] = False
  print sum(sum(testFilled))
  suboptNew=_getBestN(testFilled,filledPlex,5,reset=False)
  #print "plex", plex
  #print "subopt", subopt
  subopt = list(set(subopt + suboptNew))
#  n += 1
  
#  totalCoverage = _getTotalFills([p[-1] for p in (plex + subopt)],filled)
#  coverageComb = _getTotalFills(pIs,filled)
  if n < decs: n +=1
#  else if coverageComb == allBlocks: break
  else: 
    testFilled = _setMatrix((minRatio - n*decrementRatio), 
                          (minCount - n*decrementCount), 
                          plex + subopt)

print >>logfile, len(subopt)  
logfile.close()

#######
# print and collate results
#######

if 1==0:
  print >>sys.stdout,"#optimal"
  for p in plex:
    pattern, ratio, count, coverage, total, index = p
    print pattern, ratio, count, coverage, total, index

  print >>sys.stdout,"#suboptimal"
  for p in subopt:
    pattern, ratio, count, coverage, total, index = p
    print pattern, ratio, count, coverage, total, index





######
# order and print primers
######

indices = [p[-1] for p in plex]
patterns = [p[0] for p in plex]
indices = np.array(indices)

plexMatch = np.zeros((allBlocks,len(indices)))
pi=0
for i in indices: # fill primers in order they were picked
  plexMatch[:,pi] = matchcounts[:,i]
  pi += 1

outMatrix = np.hstack((blockPosns,plexMatch))


np.savetxt(outfile+".tab.txt", outMatrix, fmt='%s')
np.savetxt(outfile+".chrs", blockchrs, fmt='%s')
#np.savetxt(outfile+".plex", patterns, fmt='%s')
np.savetxt(outfile+".plex", plex, fmt='%s')

# print out subopt matrices
indices = [p[-1] for p in subopt]
subPatterns = [p[0] for p in subopt]
indices = np.array(indices)

subPlexMatch = np.zeros((allBlocks,len(indices)))
pi=0
for i in indices: # fill primers in order they were picked
  subPlexMatch[:,pi] = matchcounts[:,i]
  pi += 1

subOutMatrix = np.hstack((blockPosns,subPlexMatch))

np.savetxt(outfile+".sub.tab.txt", subOutMatrix, fmt='%s')
np.savetxt(outfile+".sub.chrs", blockchrs, fmt='%s')
#np.savetxt(outfile+".sub.plex", subPatterns, fmt='%s')
np.savetxt(outfile+".sub.plex", subopt, fmt='%s')


combOutMatrix = np.hstack((blockPosns,plexMatch,subPlexMatch))
combPatterns = patterns + subPatterns
combPlex = plex + subopt

np.savetxt(outfile+".comb.tab.txt", combOutMatrix, fmt='%s')
np.savetxt(outfile+".comb.chrs", blockchrs, fmt='%s')
#np.savetxt(outfile+".comb.plex", combPatterns, fmt='%s')
np.savetxt(outfile+".comb.plex", combPlex, fmt='%s')
