# created by Todd Osmundson and Thomas Roehl 2021

from Bio import SeqIO
import sys

# name of the output file in fasta format
outputFile = sys.argv[1]
# list of gene sequences and corresponding MSTRNG names
inputFile = sys.argv[2]
# list of DEGs as a text file with MSTRNG names separated by newline
degFile = sys.argv[3]

output = open(outputFile, 'a')

def seqprocess(fastafile, degListFile):
    #generate DEG list from text file
    degOpen = open(degListFile, 'r')
    degStr = degOpen.read()
    degList = degStr.splitlines()
    seqcount = 0
    fasta_sequences = SeqIO.parse(open(fastafile), 'fasta') # parse the fasta file into sequences
    for fasta in fasta_sequences:
        seqcount = seqcount + 1
        name, sequence = fasta.id, fasta.seq
        # If sequence listed in DEGs, write it to the output file:
        if name in degList:
            fastout = '>' + name + '\n' + str(sequence) + '\n'
            output.write(fastout)
    print("Total number of input sequences processed: ", seqcount)
seqprocess(inputFile, degFile)
