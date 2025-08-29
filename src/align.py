import sys
from os.path import isdir, isfile
from working.helpers import *


def main():
    sequence1 = get_fasta_seq("data\\seq1.fasta")
    sequence2 = get_fasta_seq("data\\seq2.fasta")

    ex = Align()
    print(ex.get_alignment(sequence1, sequence2))

    print("Run executed")


if __name__ == "__main__":
    main()