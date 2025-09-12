from core.utils import *
from core.needleman_wunsch import *
from core.smith_waterman import *
from flask import Flask, render_template, request, redirect
import sqlite3

connection = sqlite3.connect("data\\saved_jobs.db")
db = connection.cursor()
with open("data\\get_jobs.sql", "r") as sfile:
    sql_queries = sfile.read()
db.executescript(sql_queries)
connection.commit()
print([x for x in db.fetchall()])
connection.close()

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def index():
    return render_template("index.html",)

@app.route("/submit_form", methods=["GET", "POST"])
def submitter():
    if request.method == "POST":
        match = request.form.get("match", type=int)
        mismatch = request.form.get("mismatch", type=int)
        gap = request.form.get("gap", type=int)
        nw_seq1, nw_seq2, nw_alignment, sw_seq1, sw_seq2, sw_alignment = main(match, mismatch, gap)

    return render_template(
        "submit_form.html",
        nw_alignment=nw_alignment,
        nw_seq1=nw_seq1,
        nw_seq2=nw_seq2,
        sw_alignment=sw_alignment,
        sw_seq1=sw_seq1,
        sw_seq2=sw_seq2,
    )


def main():#match, mismatch, gap):
    sequence1_name, sequence1, is_multiseq = get_fasta_seq("data\\examples\\seq1.fasta") # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2_name, sequence2, is_multiseq = get_fasta_seq("data\\examples\\seq2.fasta") # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta
    blosum62 = get_aa_matrix("resources\\blosum62.txt")
    a = AlignNW("dna", aa_matrix=None)#, match=match, mismatch=mismatch, gap=gap)
    b = AlignSW("dna", aa_matrix=None)#, match=match, mismatch=mismatch, gap=gap)

    #print(a.get_alignment(sequence1, sequence2))
    #print(b.get_alignment(sequence1, sequence2))
    #print(sequence1_name, sequence1, sequence2_name, sequence2)

    return sequence1_name, sequence2_name, a.get_alignment(sequence1, sequence2), sequence1_name, sequence2_name, b.get_alignment(sequence1, sequence2)


if __name__ == "__main__":
    main()#match, mismatch, gap)
    #app.run(debug=True)