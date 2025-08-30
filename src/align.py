from working.helpers import *


def main():
    sequence1 = get_fasta_seq("data\\human_hbb.fasta") # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2 = get_fasta_seq("data\\puffer_hbb.fasta") # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta
    blosum62 = get_blosum_matrix("resources\\blosum62.txt")

    a = Align("protein", blosum_matrix=blosum62)
    print(a.get_alignment(sequence1, sequence2))


if __name__ == "__main__":
    main()