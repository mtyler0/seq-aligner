"""
Michael Tyler 8/22/25
"""

import numpy as np
from abc import ABC, abstractmethod


class AlignerBaseClass(ABC):

    STOP = 0
    DIAGONAL = 1
    LEFT = 2
    UP = 3

    def __init__(self, molecule, aa_matrix, match=1, mismatch=-1, gap=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.aa = aa_matrix
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

    def _create_matrix(self, seq1, seq2):
        m, n = len(seq1), len(seq2)
        scoring_matrix = np.zeros((m + 1, n + 1), dtype=int)
        path_matrix = np.zeros((m + 1, n + 1), dtype=int)

        return scoring_matrix, path_matrix, m, n

    def _traceback_vars(self):
        return "", "", "", 0, 0, 0

    def _score_cell(self, scoring_matrix, seq1, seq2, i, j):
        if self.is_nucleotide:
            diagonal_score = scoring_matrix[i - 1, j - 1] + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch)
        else:
            diagonal_score = scoring_matrix[i - 1, j - 1] + self.aa[seq1[i - 1]][seq2[j - 1]]
        vert_score = scoring_matrix[i - 1, j] + self.gap
        horizontal_score = scoring_matrix[i, j - 1] + self.gap

        return diagonal_score, vert_score, horizontal_score