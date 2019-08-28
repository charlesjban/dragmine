import re
import math
import collections
# Dragmine requires Biopython
from Bio import SeqIO
import argparse
import sys
import random
import time

# start timer
start = time.time()

# Add command line arguments
options = argparse.ArgumentParser(description='DragMine - Silk Sequence Determination')
options.add_argument("-input", help="Query sequences. Must be in standard fasta format.")
options.add_argument("-out", help="Name of output file. No suffix required")
options.add_argument("-nmer", help="Choose N-mer Length. Default is 6.", default = 6)
options.add_argument("-scramble", help="Return a profile for a randomly re-ordered sequence. 'y' or 'n' ", default = 'n')
options.add_argument("-chop", help="Choose to 'chop' sequences to obtain scores for constituent 50AA length fragments.", default = 'n')
options.add_argument("-order", help="Choose to order sequences by their evenness 'e' or network 'n' scores, or by their order in the fasta input file 'f'. Default is 'n'.", default='n')
args = options.parse_args()

# check that user has defined arguments correctly:
try:
    nmerLength = int(args.nmer)
except ValueError:
    sys.exit("nmer must be an integer")
if (args.chop != 'y' and args.chop!='n'):
    sys.exit("-chop argument must be 'y' or 'n'.")
if (args.scramble != 'y' and args.scramble != 'n'):
    sys.exit("-scramble argument must be 'y' or 'n'.")
if (args.order != 'e' and args.order != 'n' and args.order != 'f'):
    sys.exit("-order argument must be 'e'(entropy),'n' (network) or 'f'(fasta).")
# Define a class 'Sequence', which will be passed each fasta sequence (SeqIO 'fasta' object). 
# Contains basic info and profiles for normal, chopped and scrambled seqs
class Sequence:
    def __init__(self, fasta):
        self.sequence = str(fasta.seq)
        self.header = str(fasta.description).replace(",","")
        print("Analysing seq: %s"%self.header)
        self.fullProfile = Profile(self.sequence)

        if args.chop == 'y':
            print("Analysing chopped fragments.")
            self.choppedSeqs = Chopped_Seqs(self.sequence)
        if args.scramble =='y':
            print("Analysing scrambled sequence")
            self.scrambled = ScrambledSeq(self.sequence)
            self.scrambledProfile = Profile(self.scrambled.sequence)
        if args.chop == 'y' and args.scramble == 'y': 
            self.scrambledChopped = Chopped_Seqs(self.scrambled.sequence)

# Define 'Profile' class which generates each nmer, and contains the entropy and network score classes
class Profile:
    def __init__(self, sequence):
        self.sequence = sequence
        self.nmers = {}
        self.nmerList =[]
        self.length = len(self.sequence)
        for i in range(0, len(self.sequence)-nmerLength+1):
            nmer = self.sequence[i:i+nmerLength]
            if nmer not in self.nmers:
                self.nmers[nmer] = 1
            else:
                self.nmers[nmer] += 1
            self.nmerList.append(nmer)
        self.sortednmers = sorted(self.nmers.items(), reverse=True, key=lambda x: x[1])
        self.entropy = Entropy(self.length, self.nmers)
        self.network = Network(self.nmerList)

# Define the 'Network' class which calculates network score
class Network:
    def __init__(self, nmerList):
        self.edgeScore = 0
        numbnmers = len(nmerList)
        possibleEdges = numbnmers*(numbnmers-1)/2
        for i in range(len(nmerList)):
            for j in range (i + 1, len(nmerList)):
                nmer1 = nmerList[i]
                nmer2 = nmerList[j]
                self.edgeScore += match_nmers(nmer1, nmer2)
        try:
            self.normalisedEdges = math.sqrt(self.edgeScore/possibleEdges)
        # if the sequence is less than the length of the kmer, an error is generated and no score is assigned.
        except ZeroDivisionError:
            self.normalisedEdges = 0

# Define 'Entropy' class which calculates the shannons entropy (and normalised 'evenness'), after being passed the nmer profile
class Entropy:
    def __init__(self, length, nmers):
        totalnmers = length - nmerLength + 1
        H = 0
        for freq in nmers.values():
            pi = freq/totalnmers
            piLogPi = pi * math.log2(pi)
            H = H + piLogPi
        self.H = - H
        try:
            self.HNorm = self.H / math.log2(totalnmers)
        except ZeroDivisionError:
            self.HNorm = 1

# Define 'ScrambledSeq' class  which creates an amino acid string of the same length and composition as it is passed, but in a random order
class ScrambledSeq:
    def __init__(self, AAsequence):
        AAs = []
        scrambledSequence = ""
        for AA in AAsequence:
            AAs.append(AA)
        while AAs:
                randomNo = random.randint(0,len(AAs)-1)
                scrambledSequence = scrambledSequence + AAs[randomNo]
                del AAs[randomNo]
        self.sequence = scrambledSequence

# Define 'chopped  which creates contings of 50AAs with a 25AA overlap, for the length of the sequence, and returns the entropy score for each 
# note: function will always include the last 50 AAs, to ensure complete coverage, however there will be uneven coverage at the end of sequences. 
class Chopped_Seqs:
    def __init__(self, AAsequence):
        self.choppedList = []
        for i in range(0, len(AAsequence)-49, 25):
            fragment = AAsequence[i:i+50]
            fragProf = Profile(fragment)
            self.choppedList.append(fragProf)
        fragment = AAsequence[-50:]
        fragProf = Profile(fragment)
        self.choppedList.append(fragProf)

# define the function which matches the Nmers to generate the network scores.
# Arbitrary scoring (1 for complete match, 0.75 for one mismatch between kmers)
def match_nmers(nmer1, nmer2):
    count = nmerLength
    for i in range(0, nmerLength):
        if nmer1[i] != nmer2[i]:
            count -= 1
        if count == nmerLength -2:
            break
    if count == nmerLength:
        return 1
    elif count == nmerLength-1:
        return 0.75
    else:
        return 0


sequenceList = []
# for each fasta sequence create new sequence and add to list
fastaFile = args.input.strip(".csv")

for protein in SeqIO.parse(fastaFile, "fasta"):
    if len(protein.seq) >= nmerLength:
        newSeq = Sequence(protein)
        sequenceList.append(newSeq)
if args.order == 'e':
    sequenceList.sort(key=lambda x: x.fullProfile.entropy.HNorm)
elif args.order == 'n':
    sequenceList.sort(key=lambda x: x.fullProfile.entropy.HNorm, reverse=False)

numberSeqs = len(sequenceList)

outputName = args.out

print("Writing output file. Results saving to: %s%d.csv"%(outputName, nmerLength))

resultsFile = open("%s%d.csv"%(outputName, nmerLength), "w")
# print output results to file
print("Sequence, Length AAs, Entropy, Normalised Entropy, EdgeScore, Type, Highest frequency nmers", file=resultsFile)
for seq in sequenceList:
    if args.chop == 'n':
        print("%s,%s,%s,%s,%s,%s,%s"%(seq.header+" Full", seq.fullProfile.length, seq.fullProfile.entropy.H, seq.fullProfile.entropy.HNorm, seq.fullProfile.network.normalisedEdges, "RealFull", str(seq.fullProfile.sortednmers[:10]).replace(",", "")), file= resultsFile)
        if args.scramble == 'y':
            print("%s,%s,%s,%s,%s,%s,%s"%(seq.header+" Full Scrambled", seq.scrambledProfile.length, seq.scrambledProfile.entropy.H, seq.scrambledProfile.entropy.HNorm, seq.scrambledProfile.network.normalisedEdges, "ScrambledFull", str(seq.scrambledProfile.sortednmers[:10]).replace(",", "")), file= resultsFile)
    if args.chop == 'y':
        for i in range(0, len(seq.choppedSeqs.choppedList)):
            print("%s,%s,%s,%s,%s,%s,%s"%(seq.header+" Fragment %d"%(i+1), seq.choppedSeqs.choppedList[i].length, seq.choppedSeqs.choppedList[i].entropy.H, seq.choppedSeqs.choppedList[i].entropy.HNorm, seq.choppedSeqs.choppedList[i].network.normalisedEdges, "RealFrag", str(seq.choppedSeqs.choppedList[i].sortednmers[:10]).replace(",", "")), file= resultsFile)
        if args.scramble == 'y': 
            for i in range(len(seq.scrambledChopped.choppedList)):
                print("%s,%s,%s,%s,%s,%s,%s"%(seq.header+"Fragment %d Scrambled"%(i+1), seq.scrambledChopped.choppedList[i].length, seq.scrambledChopped.choppedList[i].entropy.H, seq.scrambledChopped.choppedList[i].entropy.HNorm, seq.scrambledChopped.choppedList[i].network.normalisedEdges, "ScrambledFrag", str(seq.scrambledChopped.choppedList[i].sortednmers[:10]).replace(",", "")), file= resultsFile)

# end timer and print: number of seqs processed and time taken
end = time.time()
timeTaken = end - start
print("%d sequences processed\nnmer length: %d\nTime elapsed: %s"%(numberSeqs, nmerLength, timeTaken))
