from core.aligner_base import AlignerBaseClass

class AlignSW(AlignerBaseClass):
    """
    Runs Smith-Waterman algo.
    Requires blosum62 or pam160 for proteins
    """

    def __init__(self, molecule, aa_matrix, match=1, mismatch=-1, gap=-1):
        super().__init__(molecule, aa_matrix, match=match, mismatch=mismatch, gap=gap)


    def _initialize_matrices(self, seq1, seq2):
        # Returns scoring matrix based on inputs
        return self._create_matrix(seq1, seq2)

    def populate_matrices(self, seq1, seq2):
        # Assign scores at each position in the matrix for possible movements
        scoring_matrix, path_matrix, m ,n = self._initialize_matrices(seq1, seq2)
        max_score = -1
        max_index = (-1, -1)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                diagonal_score, vert_score, horizontal_score = self._score_cell(scoring_matrix, seq1, seq2, i, j)
                score = max(0, diagonal_score, horizontal_score, vert_score)
                scoring_matrix[i, j] = score

                # Populate path matrix for backtracking based on travel score
                if score == 0:
                    path_matrix[i, j] = AlignerBaseClass.STOP
                elif score == diagonal_score:
                    path_matrix[i, j] = AlignerBaseClass.DIAGONAL
                elif score == vert_score:
                    path_matrix[i, j] = AlignerBaseClass.UP
                elif score == horizontal_score:
                    path_matrix[i, j] = AlignerBaseClass.LEFT

                # Keep track of max score and position
                if scoring_matrix[i, j] >= max_score:
                    max_score = scoring_matrix[i, j]
                    max_index = (i, j)

        return path_matrix, max_score, max_index


    def get_alignment(self, seq1, seq2):
        # Returns alignment string and score through backtracking
        # Get matrices and initialize variables to empty strings and zeros
        path_matrix, max_score, max_index = self.populate_matrices(seq1, seq2)
        top, matches, bottom, match_counter, gaps1, gaps2 = self._traceback_vars()
        (maxi, maxj) = max_index
        current_aligned1 = ""
        current_aligned2 = ""

        # Backtrack starting at bottom right of matrix and build alignment string based on path score until
        while path_matrix[maxi, maxj] != AlignerBaseClass.STOP:
            path = path_matrix[maxi, maxj]
            match_identifier = "."

            if path == AlignerBaseClass.DIAGONAL:
                current_aligned1 = seq1[maxi - 1]
                current_aligned2 = seq2[maxj - 1]
                maxi -= 1
                maxj -= 1
                if current_aligned1 == current_aligned2:
                    match_identifier = "|"
                    match_counter += 1
            elif path == AlignerBaseClass.UP:
                current_aligned2 = "-"
                gaps2 += 1
                maxi -= 1
            elif path == AlignerBaseClass.LEFT:
                current_aligned1 = "-"
                gaps1 += 1
                maxj -= 1

            # Build sequence in proper order
            top = current_aligned1 + top
            bottom = current_aligned2 + bottom
            matches = match_identifier + matches

        """
        seq1_start = seq1.index(top[0])
        seq1_end = seq1.index(top[-1])
        seq2_start = seq2.index(top[0])
        seq2_end = seq2.index(top[-1])
        leading_seq1 = seq1[:seq1_start]
        trailing_seq1 = seq1[seq1_end:]
        leading_seq2 = seq2[:seq2_start]
        trailing_seq2 = seq2[seq2_end:]
        """
        final_length = len(top)
        percent_identity = (match_counter/final_length) * 100
        gap = max(gaps1, gaps2)

        keys = ["aligned_seq1", "match_identifiers", "aligned_seq2", "percent_id", "gap_count"]
        vals = [top, matches, bottom, percent_identity, gap]

        aligned_seq = dict(zip(keys, vals))

        # Old f string for preformatted HTML insertion as plain text. Switching to data structure output to have frontend handle display of seqs
        # f"{top}\n{matches}\n{bottom}\n\nScore: {score:.0f}\nPercent Identical: {percent_identity:.1f}%\nGaps: {gap}"
        return f"{top}\n{matches}\n{bottom}\n\nScore: {max_score:.0f}\nPercent Identical: {percent_identity:.1f}%\nGaps: {gap}"#aligned_seq