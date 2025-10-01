import os
from src.core import utils

# a = needleman_wunsch.AlignNW("dna", aa_matrix=None)
# b = needleman_wunsch.AlignNW("dna", aa_matrix=None)

data_file = os.path.join(os.path.dirname(__file__), "..", "src", "data", "examples", "fw seqs r4.docx")
data_file = os.path.abspath(data_file)

def test_get_alignmentNW():
    pass

def test_get_alignmentSW():
    pass

print(utils.get_fasta_seq(data_file))

# import chardet

# with open(data_file, "rb") as f:
#     raw_data = f.read()
#     result = chardet.detect(raw_data)
#     print(result['encoding'])