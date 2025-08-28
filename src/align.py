import sys
from os.path import isdir, isfile
from working.helpers import *


def main():
    match, mismatch, gap = 1, -1, -2
    sequence1 = get_fasta_seq("data\\seq1.fasta")
    sequence2 = get_fasta_seq("data\\seq2.fasta")

    ex = Align(match, mismatch, gap)
    print(ex.get_alignment(sequence1, sequence2))
    print(ex.paths)
    
    #print(get_fasta_seq(sequence1))
    #print(get_fasta_seq(sequence2))
    print("Run executed")


if __name__ == "__main__":
    main()