from core.aligner_base import AlignerBaseClass

class AlignNW(AlignerBaseClass):
    """
    Runs Needleman-Wunsch algo
    Requires blosum62 or pam160 for proteins
    """

    def __init__(self, molecule, aa_matrix, match=1, mismatch=-1, gap=-1):
        super().__init__(molecule, aa_matrix, match=match, mismatch=mismatch, gap=gap)


    def _initialize_matrices(self, seq1, seq2):
        # Returns scoring matrix based on inputs
        scoring_matrix, path_matrix, m , n = self._create_matrix(seq1, seq2)
        
        # Initialize border rows and columns with gap penalties
        for i in range(1, m + 1):
            scoring_matrix[i, 0] = i * self.gap
            path_matrix[i, 0] = AlignerBaseClass.LEFT # Backtrack left (traveled right)
        for j in range(1, n + 1):
            scoring_matrix[0, j] = j * self.gap
            path_matrix[0, j] = AlignerBaseClass.UP # Backtrack up (traveled down)

        return scoring_matrix, path_matrix


    def populate_matrices(self, seq1, seq2):
        # Assign scores at each position in the matrix for possible movements
        scoring_matrix, path_matrix = self._initialize_matrices(seq1, seq2)
        m, n = len(seq1), len(seq2)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                diagonal_score, vert_score, horizontal_score = self._score_cell(scoring_matrix, seq1, seq2, i, j)
                score = max(diagonal_score, horizontal_score, vert_score)
                scoring_matrix[i, j] = score

                # Populate path matrix for backtracking based on travel score
                if score == diagonal_score:
                    path_matrix[i, j] = AlignerBaseClass.DIAGONAL
                elif score == vert_score:
                    path_matrix[i, j] = AlignerBaseClass.UP
                else:
                    path_matrix[i, j] = AlignerBaseClass.LEFT

        return scoring_matrix, path_matrix


    def get_alignment(self, seq1, seq2):
        # Returns alignment string and score through backtracking
        scoring_matrix, path_matrix = self.populate_matrices(seq1, seq2)
        score = scoring_matrix[-1, -1]
        top, matches, bottom, match_counter, gaps1, gaps2 = self._traceback_vars()
        i = len(seq1) - 1
        j = len(seq2) - 1

        # Backtrack starting at bottom right of matrix and build alignment string based on path score until
        while i >= 0 and j >= 0:
            path = path_matrix[i + 1, j + 1]
            current_aligned1 = seq1[i]
            current_aligned2 = seq2[j]
            match_identifier = "."

            if path == AlignerBaseClass.DIAGONAL:
                i -= 1
                j -= 1
                if current_aligned1 == current_aligned2:
                    match_identifier = "|"
                    match_counter += 1
            elif path == AlignerBaseClass.UP:
                i -= 1
                current_aligned2 = "-"
                gaps1 += 1
            elif path == AlignerBaseClass.LEFT:
                j -= 1
                current_aligned1 = "-"
                gaps2 += 1

            top = current_aligned1 + top
            bottom = current_aligned2 + bottom
            matches = match_identifier + matches

        final_length = len(top)
        if final_length < 1:
            return "No alignment possible within the given parameters"
        percent_identity = (match_counter/final_length) * 100
        gap = max(gaps1, gaps2)

        keys = ["aligned_seq1", "match_identifiers", "aligned_seq2", "percent_id", "gap_count"]
        vals = [top, matches, bottom, percent_identity, gap]

        aligned_seq = dict(zip(keys, vals))

        # Old f string for preformatted HTML insertion as plain text. Switching to data structure output to have frontend handle display of seqs
        # f"{top}\n{matches}\n{bottom}\n\nScore: {score:.0f}\nPercent Identical: {percent_identity:.1f}%\nGaps: {gap}"
        return f"{top}\n{matches}\n{bottom}\n\nScore: {score:.0f}\nPercent Identical: {percent_identity:.1f}%\nGaps: {gap}"#aligned_seq
