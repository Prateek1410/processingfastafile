#!/usr/bin/python3
import sys
from collections import Counter #to be used to add dictionaries while adding the values of common keys (relevant for finding most frequent repeat in entire file)

if len(sys.argv) != 4:
    sys.exit('Please provide all positional arguments') 

inputfastafile = sys.argv[1]
inputframenumber = int(sys.argv[2])
n = int(sys.argv[3]) #size of the repeats to be found


#(1) How many records are in the file? 

#storing all the records in a list and then finding their count
with open(inputfastafile, 'r') as f:
    records={}
    for line in f:
        if line.startswith('>') and line[1] != ' ':
            name = line.split()[0][1:]
            records[name] = ''
        else:
            line = line.strip('\n') #removing the newline character that appears at the end of a sequence
            records[name] += line
print('***ANSWER 1)***')
print(f'There are {len(records)} records in the fasta file provided')
print('\n')

#2) What are the lengths of the sequences in the file? What is the longest sequence and what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers? 

#adding all sequences into a list to find longest/shortest
seqlist = list(records.values())

longest = {}
shortest = {}
for i,j in records.items():
    if len(j) == len(max(seqlist, key=len)): #adding key and value pairs to longest dictionary if the value is equal to the longest sequence
        longest[i] = j
    if len(j) == len(min(seqlist, key=len)): #adding key and value pairs to shortest dictionary if the value is equal to the shortest sequence
        shortest[i] = j
print('***ANSWER 2)***')
print(f'Number of shortest sequence: {len(shortest.values())} and the length is {len(list(shortest.values())[0])} nt')
print(f'Number of longest sequence: {len(longest.values())} and the length is {len(list(longest.values())[0])} nt') 
print('\n')

#3) Given an input reading frame on the forward strand (1, 2, or 3) the program should be able to identify all ORFs present in each sequence of the FASTA file and: what is the length of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it?

#4) Given a length n, the program should be able to identify all repeats of length n in all sequences in the FASTA file. The program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.
#Basically, I have to find all the substrings of length 'n' occurring in any coding sequence.  

def occurrences(string, sub): #to be used for counting repeats in sequences
    
    '''Counts overlapping substrings since the built-in count method of string class doesn't consider overlap'''
    
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count
        

class codingsequence():

    def __init__(self, identifier, sequence):
        self.sequence = sequence
        self.identifier = identifier
        self.number = identifier.split('|')[-1]
        self.orf = self.findreadingframe(self.sequence, inputframenumber) #class methods to find the readingframe and the longest orf are called during the creation of objects
        self.longestorf = self.findlongestorf(self.orf)
        self.setofrepeats = {self.sequence[i:i+n] for i in range(len(self.sequence)-n+1)} #finds repeats of given length 'n'
        self.countofrepeats = {repeat:occurrences(self.sequence, repeat) for repeat in self.setofrepeats} #create a dictionary of repeat and its count/frequency
        
    def __repr__(self):
        return f"codingsequence{self.number}"
    
    def findreadingframe(self, sequence, framenumber):

        '''Given a sequence and a reading frame, this function outputs all the possible ORFs, if any, for that particular frame'''
    
        stopsequences = ['TAG','TAA','TGA']
        start,end = [],[]
        for i in range(framenumber-1, len(sequence), 3): #first frame means I'll start from the first nucleotide thus, zeroth index. 
            if sequence[i:i+3] == 'ATG':
                start.append(i)
            if sequence[i:i+3] in stopsequences:
                end.append(i)
            
        #after collecting all the start and stop codon positions I lay out the conditions for existence of an ORF. 
    
        if len(end) == 0 or len(start) == 0 or min(start) > max(end): #if there isn't any start codon OR any stop codon OR the last stop codon doesnt come after the first stop codon then there's no ORF  
            orf = None
            return orf
    
        else:
            start.sort()
            end.sort()
            orf = {}
        
            while len(start) != 0 and len(end) != 0 and min(start) < max(end): # run the loop only while ALL of these conditions are true; if any one is false, stop. 
            
                indexofstopcodon_inseq = min([x for x in end if x>min(start)]) #index of the stop codon to be used in an iteration; this codon must be the FIRST one AFTER  the end of last orf 
                orf[f'orf_{min(start)}'] = sequence[min(start):indexofstopcodon_inseq+3] #orf starting from the first unused start codon till the current stop codon; each orf has been conveniently named to include its starting position
    
                indexofstopcodon_inendlist = end.index(indexofstopcodon_inseq) #index of the stop codon WITHIN THE END LIST 
            
                #revising the start and end lists to find the next ORF, if any
                start = [x for x in start if x>end[indexofstopcodon_inendlist]] #retaining only those start codon positions that come AFTER the end of last documented ORF.
                del end[:indexofstopcodon_inendlist + 1]  #removing stop codon positions: the one that was used this iteration and the ones that come before it (if any) i.e. were skipped
        
            return orf
        
    def findlongestorf(self, orf):
            
        '''Returns the longest orf for a particular gene'''
        
        if orf!=None:
            for frame in orf:
                if len(orf[frame])==len(max(list(orf.values()), key=len)): #if a value is equal to the max length amongst the orfs (of an gene) then store that value and its key in a dictionary
                    longestorf = {frame:orf[frame]}
            return longestorf    
            
fastarecords = [codingsequence(gene, sequence) for gene,sequence in records.items()] #initialising records from fasta file as objects of class codingsequenc

all_longestorfs = [len(list(i.longestorf.values())[0]) for i in fastarecords if i.orf != None] # if orf exists for  an object/codingseq, then store its longest orf in a list
longestorfinrecords = max(all_longestorfs) #longest orf from amongst the longest orf of each object/codingseq

print('***ANSWER 3)***')

for i in fastarecords: #finding the identifier of the sequence containing the longest orf in the records
    if i.orf != None:
        if len(list(i.longestorf.values())[0]) == longestorfinrecords:
            print(f'The codingsequence having the longest orf of length {longestorfinrecords} is {i.identifier}')

print('\n')

allcountofrepeats = [Counter(fastarecords[i].countofrepeats) for i in range(len(fastarecords))]

countofrepeats_entirefile = allcountofrepeats[0]
for c in allcountofrepeats[1:]:
    countofrepeats_entirefile += c
    
countofmostfrequentrepeat = max(list(countofrepeats_entirefile.values()))

print('***ANSWER 4)***')

for repeat, count in countofrepeats_entirefile.items():
    if countofrepeats_entirefile[repeat] == countofmostfrequentrepeat:
        print(f'The most frequent repeat of length {n} is {repeat} and it occurs {countofmostfrequentrepeat} times in the entire file')
        
