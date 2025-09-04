"""
Michael Tyler 8/22/25
"""

import numpy as np
from abc import ABC, abstractmethod

def get_blosum_matrix(blosum_file):
    """
    Returns a scoring matrix from .txt file

    :param blosum_filepath: path to blosum62 matrix (for protein alignment)
    :return: nested dictionary scoring matrix
    """
    with open(blosum_file) as b:
        header = b.readline().upper().strip().replace('  ', ' ')
        species = header.split(' ')
        final_matrix = {}
        for line in b.readlines():
            row = line.upper().strip().replace('  ', ' ').split(' ')
            species2 = row[0]
            species2_score = {}
            for i in range(len(species)):
                species2_score[species[i]] = int(row[i + 1])
            final_matrix[species2] = species2_score
        return final_matrix


def get_fasta_seq(fasta_file):
    """
    Extracts a sequence from the second line of a .fasta file

    Returns:
        Sequence as string
    """

    with open(fasta_file) as fasta:
        contents = fasta.read().upper()
        sequence = contents.split("\n")[1].strip("\n")
        if "\\" in fasta_file:
            name_start = fasta_file.index("\\") + 1
            name = fasta_file[name_start:]
        else:
            name = fasta_file

        return name, sequence


class AlignerBaseClass(ABC):

    stop = 0
    diagonal = 1
    left = 2
    up = 3

    def __init__(self, molecule, blosum_matrix, match=1, mismatch=-1, gap=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.blosum = blosum_matrix
        self.is_nucleotide = molecule.lower() == "dna"


    @abstractmethod
    def _initialize_matrices(self, seq1, seq2):
        pass

    @abstractmethod
    def run_algo(self, seq1, seq2):
        pass

    @abstractmethod
    def get_alignment(self, seq1, seq2):
        pass


class AlignNW(AlignerBaseClass):
    """
    Runs Needleman-Wunsch algo
    Requires blosum62 for proteins
    """

    def __init__(self, molecule, blosum_matrix):
        super().__init__(molecule, blosum_matrix)


    def _initialize_matrices(self, seq1, seq2):
        # Returns scoring matrix based on inputs
        m, n = len(seq1), len(seq2)
        scoring_matrix = np.zeros((m + 1, n + 1))
        path_matrix = np.zeros((m + 1, n + 1))

        # Initialize border rows and columns with gap penalties
        for i in range(1, m + 1):
            scoring_matrix[i, 0] = i * self.gap
            path_matrix[i, 0] = AlignerBaseClass.left # Backtrack left (traveled right)
        for j in range(1, n + 1):
            scoring_matrix[0, j] = j * self.gap
            path_matrix[0, j] = AlignerBaseClass.up # Backtrack up (traveled down)

        return scoring_matrix, path_matrix


    def run_algo(self, seq1, seq2):
        # Assign scores at each position in the matrix for possible movements
        scoring_matrix, path_matrix = self._initialize_matrices(seq1, seq2)
        m, n = len(seq1), len(seq2)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if self.is_nucleotide:
                    diagonal_score = scoring_matrix[i - 1, j - 1] + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch)
                else:
                    diagonal_score = scoring_matrix[i - 1, j - 1] + self.blosum[seq1[i - 1]][seq2[j - 1]]
                vert_score = scoring_matrix[i - 1, j] + self.gap
                horizontal_score = scoring_matrix[i, j - 1] + self.gap
                score = max(diagonal_score, horizontal_score, vert_score)
                scoring_matrix[i, j] = score

                # Populate path matrix for backtracking based on travel score
                if score == diagonal_score:
                    path_matrix[i, j] = AlignerBaseClass.diagonal
                elif score == vert_score:
                    path_matrix[i, j] = AlignerBaseClass.up
                else:
                    path_matrix[i, j] = AlignerBaseClass.left

        return scoring_matrix, path_matrix


    def get_alignment(self, seq1, seq2):
        # Returns alignment string and score through backtracking
        scoring_matrix, path_matrix = self.run_algo(seq1, seq2)
        score = scoring_matrix[-1, -1]
        top = "" # Reference strand w/ gaps
        matches = "" # Visual confirmation of match
        bottom = "" # Query strand w/ gaps
        i = len(seq1) - 1
        j = len(seq2) - 1
        match_counter = 0
        gaps = 0
        gaps2 = 0

        # Backtrack starting at bottom right of matrix and build alignment string based on path score until
        while (i, j) != (-1, -1):
            path = path_matrix[i + 1, j + 1]
            current_aligned1 = seq1[i]
            current_aligned2 = seq2[j]
            match_identifier = " "

            if path == AlignerBaseClass.diagonal:
                i -= 1
                j -= 1
                if current_aligned1 == current_aligned2:
                    match_identifier = "|"
                    match_counter += 1
            elif path == AlignerBaseClass.up:
                i -= 1
                current_aligned2 = "-"
                gaps += 1
            elif path == AlignerBaseClass.left:
                j -= 1
                current_aligned1 = "-"
                gaps2 += 1

            top = current_aligned1 + top
            bottom = current_aligned2 + bottom
            matches = match_identifier + matches
            gap = max(gaps, gaps2)

        final_length = len(top)
        percent_identity = (match_counter/final_length) * 100

        return f"{top}\n{matches}\n{bottom}\n\nScore: {score:.0f}\nPercent Identical: {percent_identity:.1f}%\nGaps: {gap}"


class AlignSW(AlignerBaseClass):
    """
    Runs Smith-Waterman algo.
    Requires blosum62 for proteins
    """

    def __init__(self, molecule, blosum_matrix):
        super().__init__(molecule, blosum_matrix)


    def _initialize_matrices(self, seq1, seq2):
        # Returns scoring matrix based on inputs
        m, n = len(seq1), len(seq2)
        scoring_matrix = np.zeros((m + 1, n + 1), dtype=int)
        path_matrix = np.zeros((m + 1, n + 1), dtype=int)

        return scoring_matrix, path_matrix


    def run_algo(self, seq1, seq2):
        # Assign scores at each position in the matrix for possible movements
        scoring_matrix, path_matrix = self._initialize_matrices(seq1, seq2)
        m, n = len(seq1), len(seq2)
        max_score = -1
        max_index = (-1, -1)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if self.is_nucleotide:
                    diagonal_score = scoring_matrix[i - 1, j - 1] + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch)
                else:
                    diagonal_score = scoring_matrix[i - 1, j - 1] + self.blosum[seq1[i - 1]][seq2[j - 1]]
                vert_score = scoring_matrix[i - 1, j] + self.gap
                horizontal_score = scoring_matrix[i, j - 1] + self.gap
                score = max(0, diagonal_score, horizontal_score, vert_score)
                scoring_matrix[i, j] = score

                # Populate path matrix for backtracking based on travel score
                if score == 0:
                    path_matrix[i, j] = AlignerBaseClass.stop
                elif score == diagonal_score:
                    path_matrix[i, j] = AlignerBaseClass.diagonal
                elif score == vert_score:
                    path_matrix[i, j] = AlignerBaseClass.up
                elif score == horizontal_score:
                    path_matrix[i, j] = AlignerBaseClass.left

                # Keep track of max score and position
                if scoring_matrix[i, j] >= max_score:
                    max_score = scoring_matrix[i, j]
                    max_index = (i, j)

        return path_matrix, max_score, max_index


    def get_alignment(self, seq1, seq2):
        # Returns alignment string and score through backtracking
        path_matrix, max_score, max_index = self.run_algo(seq1, seq2)
        top = "" # Reference strand w/ gaps
        matches = "" # Visual confirmation of match
        bottom = "" # Query strand w/ gaps
        match_counter = 0
        gaps = 0
        gaps2 = 0
        (maxi, maxj) = max_index
        current_aligned1 = ""
        current_aligned2 = ""

        # Backtrack starting at bottom right of matrix and build alignment string based on path score until
        while path_matrix[maxi, maxj] != AlignerBaseClass.stop:
            path = path_matrix[maxi, maxj]
            match_identifier = " "

            if path == AlignerBaseClass.diagonal:
                current_aligned1 = seq1[maxi - 1]
                current_aligned2 = seq2[maxj - 1]
                maxi -= 1
                maxj -= 1
                match_identifier = "|"
                match_counter += 1
            elif path == AlignerBaseClass.up:
                current_aligned2 = "-"
                gaps2 += 1
                maxi -= 1
            elif path == AlignerBaseClass.left:
                current_aligned1 = "-"
                gaps += 1
                maxj -= 1
            
            # Build sequence in proper order
            top = current_aligned1 + top
            bottom = current_aligned2 + bottom
            matches = match_identifier + matches

        seq1_start = seq1.index(top[0])
        seq1_end = seq1.index(top[-1])
        seq2_start = seq2.index(top[0])
        seq2_end = seq2.index(top[-1])
        leading_seq1 = seq1[:seq1_start]
        trailing_seq1 = seq1[seq1_end:]
        leading_seq2 = seq2[:seq2_start]
        trailing_seq2 = seq2[seq2_end:]
        final_length = len(top)
        percent_identity = (match_counter/final_length) * 100

        return f"{top}\n{matches}\n{bottom}\n\nScore: {max_score:.0f}\nPercent Identical: {percent_identity:.1f}% ({match_counter}/{final_length})"