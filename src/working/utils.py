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
        contents = fasta.read().upper()
        sequence = contents.split("\n")[1].strip("\n")
        if "\\" in fasta_file:
            name_start = fasta_file.index("\\") + 1
            name = fasta_file[name_start:]
        else:
            name = fasta_file

        return name, sequence
