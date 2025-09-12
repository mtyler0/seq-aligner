def get_aa_matrix(matrix_file):
    """
    Returns a scoring matrix from .txt file

    :param aa_filepath: path to blosum62 or pam160 matrix (for protein alignment)
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
        seqs = []
        current_sequence = ""
        is_multiseq = False
        for i, line in enumerate(contents):
            if line[0] == ">":
                current_sequence = contents[i + 1].strip()
                seqs.append(current_sequence)

        # Fix!!! Disallow files with >2 seqs?
        if "\\" in fasta_file:
            name_start = fasta_file.index("\\") + 1
            name_end = fasta_file.index(".")
            name = fasta_file[name_start:name_end]
        else:
            name = fasta_file

        if len(seqs) == 1:
            return name, seqs[0], is_multiseq
        elif len(seqs) == 0:
            return "Empty File"
        else:
            is_multiseq = True
            return name, seqs, is_multiseq