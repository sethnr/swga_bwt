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
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


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
forwardWeight = args.forwardWeight
outfile = args.out


k = args.plexSize
minRatio = args.minRatio # 10
minCount = args.minCount # 10


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
ratios = matchcountsHDF5['ratios']
emptyBlocks = np.array([True] * allBlocks)

#filled = np.greater(matchcounts[:,:20],[0])
filled = np.greater(matchcounts[:,:],[0])
sdcols = np.std(matchcounts,axis=0)
print >>sys.stderr, filled.shape, ratios.shape, sdcols.shape

passRatio = np.greater_equal(ratios, [minRatio])
failRatio = np.less(ratios, [minRatio])
filled[:,failRatio] = False #remove those which are too low ratio

patternCounts = np.array(matchcounts).sum(axis=0)
passCount = np.greater_equal(patternCounts, [minCount])
failCount = np.less(patternCounts, [minCount])
filled[:,failCount] = False #remove those which are too low count


#filledTot = filled.sum(axis=0)

filledPlex = np.zeros((allBlocks,),dtype=bool)

#internal method to reset matrix - assumes filled/matchcounts/etc already there
def _setMatrix(passRatio, passCount, plex):
  print >>sys.stderr, "resetting filled matrix"
  testFilled = np.copy(filled)

  passRatio = np.greater_equal(ratios, [minRatio])
  failRatio = np.less(ratios, [minRatio])
  testFilled[:,failRatio] = False #remove those which are too low ratio

  patternCounts = np.array(matchcounts).sum(axis=0)
  passCount = np.greater_equal(patternCounts, [minCount])
  failCount = np.less(patternCounts, [minCount])
  testFilled[:,failCount] = False #remove those which are too low count

  filledPlex = np.zeros((allBlocks,),dtype=bool)
  indices = [p[-1] for p in plex]
#    print >> sys.stderr, indices
  if (len(indices) > 0):
    testFilled[:,np.array(indices)] = 0
  return testFilled

testFilled = _setMatrix(passRatio, passCount, [])

plex=[]
pi = 0
while pi < k:
  filledTot = testFilled.sum(axis=0)
  nextBestI = filledTot.argsort()[-1]
#  filledPlex
#  filled[:,nextBestI]
  filledPlex[testFilled[:,nextBestI]] = True
  print pi, sum(filledPlex), sum(sum(testFilled)), patternCounts[nextBestI], patterns[nextBestI]

  testFilled[filledPlex,:] = False
  plex.append(
    (patterns[nextBestI],ratios[nextBestI],sum(filled[:,nextBestI]),patternCounts[nextBestI],nextBestI))
  pi += 1
  if sum(sum(testFilled))==0:
#    print >>sys.stderr, "resetting filled matrix"
#    testFilled = np.copy(filled)
#    filledPlex = np.zeros((allBlocks,),dtype=bool)
#    indices = [p[-1] for p in plex]
#    print >> sys.stderr, indices
#    testFilled[:,np.array(indices)] = 0
    testFilled = _setMatrix(passRatio, passCount, plex)

for p in plex:
  pattern, ratio, coverage, total, index = p
  print pattern, ratio, coverage, total, index


indices = [p[-1] for p in plex]
patterns = [p[0] for p in plex]
#print >>sys.stderr, indices
#reindex = np.argsort(indices)
print >>sys.stderr, allBlocks
#indices.sort()

#patterns = [patterns[i] for i in reindex]
indices = np.array(indices)
#plexMatch = matchcounts[:,indices]
plexMatch = np.zeros((allBlocks,len(indices)))
pi=0
for i in indices:
  print matchcounts[:,i]
  plexMatch[:,pi] = matchcounts[:,i]
  pi += 1
#plexMatch = plexMatch.astype(int)
#outfile = outfile+".k"+str(k)+".r"+str(minRatio)+".tab.txt"
#posMatrix = np.hstack((blockchrs[:,1],blockPosns))
#print posMatrix.shape, plexMatch.shape
outMatrix = np.hstack((blockPosns,plexMatch))
if minRatio % 1 == 0: minRatio = int(minRatio)

outfileVars = outfile+"_k"+str(k)+"_r"+str(minRatio)+"_c"+str(minCount)
np.savetxt(outfileVars+".tab.txt", outMatrix, fmt='%s')
np.savetxt(outfileVars+".chrs", blockchrs, fmt='%s')
np.savetxt(outfileVars+".plex", patterns, fmt='%s')
