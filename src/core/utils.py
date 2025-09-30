def get_aa_matrix(matrix_file):
    """
    Returns a scoring matrix from .txt file

    :param matrix_file: path to blosum62 or pam160 matrix (for protein alignment)
    :return: nested dictionary scoring matrix
    """
    with open(matrix_file) as b:
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
        contents = fasta.readlines()
        extension = fasta_file.rsplit(".", 1)[1].lower()
        if extension == "fasta":
            seq = contents[1:]
            seq = "".join(seq).strip()
        elif extension in ("docx", "txt"):
            seq = "".join(contents)

        # Get name from file
        if "\\" in fasta_file:
            name_start = fasta_file.index("\\") + 1
            name_end = fasta_file.index(".")
            name = fasta_file[name_start:name_end]
        else:
            name = fasta_file

        return name, seq


def format_alignment(seq1, matches, seq2, line_length=60, start1=1, start2=1):
        """
        Formats alignment in an NCBI BLAST-like style:
        seq1 and seq2 are aligned strings (with gaps),
        matches is the '|' or '.' line.
        """
        output_lines = []
        i1, i2 = start1, start2  # Current positions in original sequences

        for i in range(0, len(seq1), line_length):
            s1_chunk = seq1[i:i+line_length]
            s2_chunk = seq2[i:i+line_length]
            m_chunk = matches[i:i+line_length]

            # Figure out start and end positions without counting gaps
            start_pos1 = i1
            start_pos2 = i2
            end_pos1 = i1 + sum(c != "-" for c in s1_chunk) - 1
            end_pos2 = i2 + sum(c != "-" for c in s2_chunk) - 1

            # Build each section
            top_line = f"Query  {start_pos1:<4} {s1_chunk}    {end_pos1+1 if end_pos1>=start_pos1 else start_pos1}"
            mid_line = f"       {'':<4} {m_chunk}"
            bot_line = f"Sbjct  {start_pos2:<4} {s2_chunk}    {end_pos2+1 if end_pos2>=start_pos2 else start_pos2}"

            output_lines.extend([top_line, mid_line, bot_line, ""])

            # Update sequence counters
            i1 = end_pos1 + 1
            i2 = end_pos2 + 1

        return "\n".join(output_lines)