from src.core import needleman_wunsch as nw
from src.core import smith_waterman as sw
from src.core.utils import *


def test_alignmentNW():
    n = nw.AlignNW("DNA", None, match=1, mismatch=-1, gap=-2)
    result = n.get_alignment("GATTCA","GACTGU")
    score = result[1]
    assert score == 0
    
    
def test_fasta_parsing():
    seq_name, seq = get_fasta_seq("seq1.fasta")
    assert seq_name == ">Sequence1:ABCD1"
    assert seq == "AAAAATTTCTCGGCGTCCGCGCTAA"

    
if __name__ == "__main__":
    test_alignmentNW()