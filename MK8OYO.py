'''Author: Matt Kachlik.
This program uses shotgun sequencing and utilizes both the TSP algorithm to find 
the contiguous sequence
02/24/2015
'''
import random
from random import randint
import math
import collections
from collections import defaultdict
'''
preSuffixList
makes a list for which sequences are prefixes and which are suffixes when they were aligned based of
the greedy algorithm to find the highest
'''
def preSuffixList(overLapDict,numFrag):
	prefixList = []
	suffixList = []
	high = 0
	pre = ""
	suff = ""
	for numFrag in range(0,numFrag-1):
		high = 0
		pre = ""
		suff = ""
		for k,v in overLapDict.items():
			for k2,v2 in v.items():
				if(overLapDict[k][k2] > high):
					high = overLapDict[k][k2]
					pre = k
					suff = k2
 		prefixList.append(pre)
		suffixList.append(suff)


		#removing all incoming for the suffix
		for k,v in overLapDict.items():
			for k2,v2 in v.items():
				if(k2 == suff):
					del overLapDict[k][k2]
		for k,v in overLapDict.items():
			if(k == pre):
				del overLapDict[k]
	#finding the instance(s) where there is a prefix but its not in the suffix
	start = set(prefixList) - set(suffixList)
	numStart = len(start)


	return prefixList,suffixList,start

'''
coverageMet
checks if the coverage of all the sequences are met by checking
against the fold
'''
def coverageMet(coverage, fold):
	i = 0
	met = True
	while (i < len(coverage) and met == True):
		if(coverage[i] < fold):
			met = False
		i += 1
	return met

'''
overlap
determines if the fragments of the lists that were sent
have any overlap then returns a dictionary containing all the overLaps
and how much they overLap
'''
def overlap (numFrag,frags):
	overLapDict = collections.defaultdict(lambda:defaultdict(int))
	frag1 = ""
	frag2 = ""
	overLap = 0
	#starting from the highest number of overlap to lowest possible overlap 
	for i in range(0,numFrag):
		for j in range(i+1, numFrag):
			frag1 = frags[i]
			frag2 = frags[j]
			f1Len = len(frag1)
			f2Len = len(frag2)
			#comparing suffix of frag1 to prefix of frag2
			overLap,overLapDict,frag1,frag2 = findOverLap(frag1,frag2)
			if(overLap != 0):
				if(frag1[f1Len-overLap: f1Len] == frag2[0:overLap]):
				 	overLapDict[frags.index(frag1)][frags.index(frag2)] = overLap
				else:
				 	overLapDict[frags.index(frag2)][frags.index(frag1)] = overLap
			else:
				print frag1+ " 1"
				print frag2+ " 2"
		
	return overLapDict	
'''
findOverLap
finds and returns the overlap for two frags
'''
def findOverLap(frag1,frag2):

	f1Len = len(frag1)
	f2Len = len(frag2)
	minLen = min(f1Len, f2Len)
	overlap = 0

	k = (minLen - 1)

	if(frag1 not in frag2 and frag2 not in frag1 and frag1 != frag2):
		while (k >= 1 and overlap == 0):
			if (frag1[f1Len-k: f1Len] == frag2[0:k]):
				overlap = k
				overLapDict[frags.index(frag1)][frags.index(frag2)] = overlap
				#print frag1, frag2, overLapDict, overlap			
			elif (frag2[f2Len-k:f2Len]==frag1[0:k]):
				overlap = k
				overLapDict[frags.index(frag2)][frags.index(frag1)] =  overlap
				#print frag1, frag2, overLapDict, overlap
			k -= 1
	return overlap,overLapDict,frag1,frag2

#input
overLapDict = collections.defaultdict(lambda:defaultdict(int))
#input
if1 = open(raw_input('What is the name of the fasta .txt file of the sequence?\n'))
seq = ""
for line in if1:
	line = line.replace('\n', '')
	line = line.replace('\r','')
	if(line[0] != ">"):
		seq += line

minimum = int(raw_input("Whats the min of sequence?\n"))
maximum = int(raw_input("Whats the max of sequence?\n"))
cFold = int(raw_input("What is the coverage you would like?\n"))
coverage = collections.defaultdict(int)
for i in range(0,len(seq)-1):
	coverage[i] = 0
#step #1: Generate a set fragments from the input seq
numFrag = 0
frags = []
# looping through until our coverage is met
while (coverageMet(coverage, cFold) == False): 
	randL = randint(minimum,maximum)
	randS = randint(0,len(seq)-randL)
	newFrag = seq[randS:randS+randL]

	frags.append(seq[randS:randS+randL])
	numFrag += 1
	#update coverage
	for i in range(randS, randS+randL -1):
		coverage[i] += 1


	#step 2: determine overlap
overLapDict = overlap(numFrag,frags)
prefixList,suffixList,start = preSuffixList(overLapDict,numFrag)

#iterating through start points(when is multiple) to find how the sequence lines up using the TSP algorithm 
seq = ""
print start
if(len(start) > 1):
	print "no optimal solution found printing possible contiguous:"
else:
	print "optimal solution:"
for i in start:
	print "sequence for " + str(i) + " as start point"
	suff = suffixList[0]
	pre = i 
	suff = suffixList[prefixList.index(i)]
	inSuffixList = True
	#using a boolean to see if the new suffix exists or if its the end
	while inSuffixList:
		overLap, overLapDict, frag1,frag2  = findOverLap(frags[pre],frags[suff])
		frag1Len = len(frag1)
		frag2Len = len(frag2)
		if(frag1[frag1Len-overLap: frag1Len] == frag2[0:overLap]):
			seq += (frag1[0: frag1Len - overLap] + frag2)
		elif frag2[frag2Len-overLap:frag2Len]== frag1[0:overLap]:
			seq += (suff[0:frag2Len-frag1Lap] + frag1)
		pre = suff
		#once you find that you are unable to make the suff since its not in the list leave the loop
		try:
	  		suff = suffixList[prefixList.index(pre)]
		except: 
	  		inSuffixList = False
	print seq
			

	
