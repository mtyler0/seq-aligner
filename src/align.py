from working.helpers import *
from flask import Flask, render_template


app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html", 
                           nw_alignment=main()[2], 
                           nw_seq1=main()[0], 
                           nw_seq2=main()[1],
                           sw_alignment=main()[5],
                           sw_seq1=main()[3],
                           sw_seq2=main()[4],
                           )


def main():
    sequence1_name, sequence1 = get_fasta_seq("data\\seq1.fasta") # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2_name, sequence2 = get_fasta_seq("data\\seq2.fasta") # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta
    blosum62 = get_blosum_matrix("resources\\blosum62.txt")
    a = AlignNW("dna", blosum_matrix=blosum62)
    b = AlignSW("dna", blosum_matrix=blosum62)

    #print(a.get_alignment(sequence1, sequence2))
    print(b.get_alignment(sequence1, sequence2))

    #return sequence1_name, sequence2_name, a.get_alignment(sequence1, sequence2), sequence1_name, sequence2_name, b.get_alignment(sequence1, sequence2)
    #return sequence1_name, sequence2_name, b.get_alignment(sequence1, sequence2)

if __name__ == "__main__":
    main()
    #app.run(debug=True)