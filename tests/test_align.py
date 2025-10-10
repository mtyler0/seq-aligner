from src.core import needleman_wunsch as nw
from src.core import smith_waterman as sw
from src.core import utils


def test_alignmentNW():
    n = nw.AlignNW("DNA", None, match=1, mismatch=-1, gap=-2)
    result = n.get_alignment("GATTCA","GACTGU")
    # assert result[1] > 0
    # assert "G" in result[0][0]
    print(result)


if __name__ == "__main__":
    test_alignmentNW()