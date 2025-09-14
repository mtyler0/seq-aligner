from core.utils import *
from core.needleman_wunsch import *
from core.smith_waterman import *
from flask import Flask, render_template, request, redirect, g
import sqlite3


app = Flask(__name__)
DATABASE = "data\\saved_jobs.db"

def get_db():
    db = getattr(g, "_database", None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
        db.row_factory = sqlite3.Row
        return db


@app.teardown_appcontext
def close_conneciton(exception):
    db = getattr(g, "_database", None)
    if db is not None:
        db.close()


@app.route("/", methods=["GET", "POST"])
def index():
    return render_template("index.html")


@app.route("/submit_form", methods=["GET", "POST"])
def submit():
    if request.method == "POST":
        match = request.form.get("match", type=int)
        mismatch = request.form.get("mismatch", type=int)
        gap = request.form.get("gap", type=int)
        nw_seq1, nw_seq2, nw_alignment, sw_seq1, sw_seq2, sw_alignment = params = main(match, mismatch, gap)

        db = get_db()
        cursor = db.cursor() # type: ignore
        try:
            cursor.execute(
                """
                INSERT INTO jobs (
                nw_seq1_name, 
                nw_seq2_name, 
                nw_result, 
                sw_seq1_name, 
                sw_seq2_name, 
                sw_result
                ) 
                VALUES (?, ?, ?, ?, ?, ?)
                """, params
                )
            cursor.execute("INSERT INTO jobs (match_score, mismatch, gap) VALUES (?, ?, ?)", (match, mismatch, gap))
            db.commit() # type: ignore
            job_id = cursor.lastrowid
            print(job_id)

        except Exception as e:
            db.rollback() # type: ignore
            print("DB insert failed:", e)

    return render_template(
        "submit_form.html",
        nw_alignment=nw_alignment,
        nw_seq1=nw_seq1,
        nw_seq2=nw_seq2,
        sw_alignment=sw_alignment,
        sw_seq1=sw_seq1,
        sw_seq2=sw_seq2,
    )

@app.route("/jobs")
def get_jobs():
    db = get_db()
    cursor = db.cursor() # type: ignore
    cursor.execute("SELECT rowid, * FROM jobs ORDER BY current_time DESC")
    results = cursor.fetchall()
    return render_template("jobs.html",results=results)


def main(match, mismatch, gap):
    sequence1_name, sequence1, is_multiseq = get_fasta_seq("data\\examples\\seq1.fasta") # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2_name, sequence2, is_multiseq = get_fasta_seq("data\\examples\\seq2.fasta") # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta
    blosum62 = get_aa_matrix("resources\\blosum62.txt")
    a = AlignNW("dna", aa_matrix=None, match=match, mismatch=mismatch, gap=gap)
    b = AlignSW("dna", aa_matrix=None, match=match, mismatch=mismatch, gap=gap)

    #print(a.get_alignment(sequence1, sequence2))
    #print(b.get_alignment(sequence1, sequence2))
    #print(sequence1_name, sequence1, sequence2_name, sequence2)

    return sequence1_name, sequence2_name, a.get_alignment(sequence1, sequence2), sequence1_name, sequence2_name, b.get_alignment(sequence1, sequence2)


if __name__ == "__main__":
    #main()#match, mismatch, gap)
    app.run(debug=True)