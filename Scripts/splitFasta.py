from os import sys
import math


#Usage:
#python <this script> <fasta file> <split length>

from Bio import SeqIO


splitLen = int(sys.argv[2])

fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
for fasta in fasta_sequences:
  name, sequence = fasta.id, str(fasta.seq)
  for i in range(0, int(math.ceil(1.0 * len(sequence)/splitLen))):
    new_name = ">" + name + "_" + str(i)
    new_sequence = sequence[i*splitLen : (i+1)*splitLen]
    print new_name
    print new_sequence
