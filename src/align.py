from working.utils import *
from working.needleman_wunsch import *
from working.smith_waterman import *
from flask import Flask, render_template, request


app = Flask(__name__)

@app.route("/", methods=["GET, POST"])
def index():
    if request.method == "POST":
        match = request.form["match"]
        mismatch = request.form["mismatch"]
        gap = request.form["gap"]

    nw_seq1, nw_seq2, nw_alignment, sw_seq1, sw_seq2, sw_alignment = main()

    return render_template(
        "index.html",
        nw_alignment=nw_alignment,
        nw_seq1=nw_seq1,
        nw_seq2=nw_seq2,
        sw_alignment=sw_alignment,
        sw_seq1=sw_seq1,
        sw_seq2=sw_seq2,
        m1=match,
        m2=mismatch,
        g=gap
    )


def main():
    sequence1_name, sequence1 = get_fasta_seq("data\\seq1.fasta") # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2_name, sequence2 = get_fasta_seq("data\\seq2.fasta") # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta
    blosum62 = get_aa_matrix("resources\\blosum62.txt")
    a = AlignNW("dna", aa_matrix=None)
    b = AlignSW("dna", aa_matrix=None)

    print(a.get_alignment(sequence1, sequence2))
    print(b.get_alignment(sequence1, sequence2))

    return sequence1_name, sequence2_name, a.get_alignment(sequence1, sequence2), sequence1_name, sequence2_name, b.get_alignment(sequence1, sequence2)
    #return sequence1_name, sequence2_name, b.get_alignment(sequence1, sequence2)

if __name__ == "__main__":
    #main()
    app.run(debug=True)