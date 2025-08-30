from working.helpers import *
from flask import Flask, render_template


app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html", alignment=main()[2], seq1=main()[0], seq2=main()[1])


def main():
    sequence1_name, sequence1 = get_fasta_seq("data\\seq1.fasta") # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2_name, sequence2 = get_fasta_seq("data\\seq2.fasta") # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta
    blosum62 = get_blosum_matrix("resources\\blosum62.txt")
    a = Align("dna", blosum_matrix=blosum62)

    return sequence1_name, sequence2_name, a.get_alignment(sequence1, sequence2)


if __name__ == "__main__":
    app.run(debug=True)