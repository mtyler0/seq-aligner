"""

8/22/25
Michael Tyler

"""

import numpy as np

def get_fasta_seq(fasta_file):
    """
    Extracts a sequence from the second line of a .fasta file

    Returns:
        Sequence as string
    """

    with open(fasta_file) as fasta:
        contents = fasta.read().upper()
        sequence = contents.split("\n")[1].strip("\n")
        return sequence

class Align:
    """
    Runs Needleman-Wunsch algo
    Requires scoring matrix
    """

    def __init__(self, match=1, mismatch=-1, gap=-2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap


    def get_alignment(self, seq1, seq2):
        # Returns alignment string and score through backtracking
        scoring_matrix, path_matrix = self.needleman_algo(seq1, seq2)
        score = scoring_matrix[-1][-1]
        top = "" # Reference strand w/ gaps
        matches = "" # Visual confirmation of match
        bottom = "" # Query strand w/ gaps
        i = len(seq1) - 1
        j = len(seq2) - 1

        # Backtrack starting at bottom right of matrix and build alignment string based on path score until
        while (i, j) != (-1, -1):
            path = path_matrix[i+1][j+1]
            letter1 = seq1[i]
            letter2 = seq2[j]
            match_identifier = " "

            if path == 0:
                i -= 1
                j -= 1
                if letter1 == letter2:
                    match = "|"
            elif path == 2:
                i -= 1
                letter2 = "-"
            else:
                j -= 1
                letter1 = "-"

            top = letter1 + top
            bottom = letter2 + bottom
            matches = match_identifier + matches

        return f"{top}\n{matches}\n{bottom}\nScore: {score}"


    def _initialize_matrices(self, seq1, seq2):
        # Returns scoring matrix based on inputs
        m, n = len(seq1), len(seq2)
        scoring_matrix = np.zeros((m+1, n+1))
        path_matrix = np.zeros((m+1, n+1))

        # Initialize border rows and columns with gap penalties
        for i in range(1, m+1):
            scoring_matrix[i][0] = i * self.gap
            path_matrix[i][0] = 1 # Backtrack left (traveled right)
        for j in range(1, n+1):
            scoring_matrix[0][j] = j * self.gap
            path_matrix[0][j] = 2 # Backtrack up (traveled down)
        return scoring_matrix, path_matrix


    def needleman_algo(self, seq1, seq2):
        # Assign scores at each position in the matrix for possible movements
        scoring_matrix, path_matrix = self._initialize_matrices(seq1, seq2)
        m, n = len(seq1), len(seq2)

        for i in range(1, m+1):
            for j in range(1, n+1):
                diagonal_score = scoring_matrix[i-1][j-1] + (self.match if seq1[i-1] == seq2[j-1] else self.mismatch)
                down_score = scoring_matrix[i-1][j] + self.gap
                right_score = scoring_matrix[i][j-1] + self.gap
                score = max(diagonal_score, right_score, down_score)
                scoring_matrix[i][j] = score

                # Populate path matrix for backtracking based on travel score
                if score == diagonal_score:
                    path_matrix[i][j] = 0
                elif score == down_score:
                    path_matrix[i][j] = 2
                else:
                    path_matrix[i][j] = 1

        self.path_matrix = path_matrix
        self.scoring_matrix = scoring_matrix

        return scoring_matrix, path_matrix