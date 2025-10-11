from src.core import needleman_wunsch as nw
from src.core import smith_waterman as sw
from src.core.utils import *


def test_alignmentNW():
    n = nw.AlignNW("DNA", None, match=1, mismatch=-1, gap=-2)
    s = sw.AlignSW("DNA", None, match=1, mismatch=-1, gap=-2)
    result_n = n.get_alignment("GATTCA","GACTGU")
    result_s = s.get_alignment("GATTCA","GACTGU")
    score_n = result_n[1]
    score_s = result_s[1]
    assert score_n == 0
    assert score_s == 2
    
    
def test_fasta_parsing():
    seq_name, seq = get_fasta_seq("seq1.fasta")
    assert seq_name == ">Sequence1:ABCD1"
    assert seq == "AAAAATTTCTCGGCGTCCGCGCTAA"

    
if __name__ == "__main__":
    test_alignmentNW()