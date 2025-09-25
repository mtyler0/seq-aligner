from core.utils import *
from core.needleman_wunsch import *
from core.smith_waterman import *
from flask import Flask, render_template, request, redirect, g, flash, url_for
import sqlite3
import os
from werkzeug.utils import secure_filename


app = Flask(__name__)
app.secret_key = "2nf38-f2n398"
DATABASE = "data\\saved_jobs.db"
UPLOAD_FOLDER = "data\\uploaded_files"
ALLOWED_EXTENSIONS = {"fasta", "txt", "docx"}
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER


def get_db():
    db = getattr(g, "_database", None)
    if db is None:
        db = g._database = sqlite3.connect(DATABASE)
        db.row_factory = sqlite3.Row
        return db


def allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def save_file(file):
    if not file.filename:
        flash("No selected file")
        return redirect("/")
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        path = os.path.join(app.config["UPLOAD_FOLDER"], filename)
        file.save(path)
        return path


@app.teardown_appcontext
def close_conneciton(exception):
    db = getattr(g, "_database", None)
    if db is not None:
        db.close()


@app.route("/", methods=["GET", "POST"])
def index():
    return render_template("index.html")


# Send data to algorithm, store in db, then redirect to submit page
@app.route("/submit_form", methods=["POST"])
def post_form():
    text_input: str = request.form["SUBJECT"].strip() and request.form["QUERY"].strip()
    input: bool = len(text_input) > 0 or \
        (request.files.get("seq1file").filename != "" and \
            request.files.get("seq2file").filename != "")
    match = request.form.get("match", type=int)
    mismatch = request.form.get("mismatch", type=int)
    gap = request.form.get("gap", type=int)
    molecule = request.form.get("type")
    matrix = request.form.get("protein-matrix")
    
    if not input:
        flash("Error: No file/text input")
        return redirect("/")
    elif text_input:
        seq1_text = request.form["SUBJECT"]
        seq2_text = request.form["QUERY"]
        try:
            params = main(match, mismatch, gap, seq1_text, seq2_text, molecule, is_text=True)
        except Exception as e:
            flash(f"ERROR: {e}")
            return redirect("/")
    else:
        seq1_file = request.files["seq1file"]
        seq2_file = request.files["seq2file"]
        seq1_path = save_file(seq1_file)
        seq2_path = save_file(seq2_file)
        params = main(match, mismatch, gap, seq1_path, seq2_path, molecule)

    db = get_db()
    cursor = db.cursor() # type: ignore
    try:
        cursor.execute(
            """
            INSERT INTO jobs (
            match_score,
            mismatch,
            gap,
            nw_seq1_name, 
            nw_seq2_name, 
            nw_result, 
            sw_seq1_name, 
            sw_seq2_name, 
            sw_result,
            score,
            percent_id,
            gaps,
            score2,
            percent_id2,
            gaps2
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (match, mismatch, gap, *params)
            )
        db.commit() # type: ignore
        job_id = cursor.lastrowid
        # flash("Submitting alignment job...")
        return redirect(url_for("submit", job_id=job_id))
    
    except Exception as e:
        db.rollback() # type: ignore
        print("DB insert failed:", e)
        return redirect(url_for("index"))


# Get alignment data for the current job ID
@app.route("/submit")
def submit():
    job_id = request.args.get("job_id")
    db = get_db()
    cursor = db.cursor() # type: ignore
    cursor.execute("SELECT * FROM jobs WHERE rowid = ?", (job_id,))
    job = cursor.fetchone()

    if job is None:
        flash("Job not found.", "error")
        return redirect(url_for("index"))

    return render_template("submit_form.html", 
        nw_alignment=job["nw_result"],
        nw_seq1=job["nw_seq1_name"],
        nw_seq2=job["nw_seq2_name"],
        nw_score=job["score"],
        nw_percent_id=job["percent_id"],
        nw_gaps=job["gaps"],
        sw_alignment=job["sw_result"],
        sw_score=job["score2"],
        sw_percent_id=job["percent_id2"],
        sw_gaps=job["gaps2"],
        sw_seq1=job["sw_seq1_name"],
        sw_seq2=job["sw_seq2_name"])


# Retrieve all past data
@app.route("/jobs")
def get_jobs():
    db = get_db()
    cursor = db.cursor() # type: ignore
    cursor.execute("SELECT rowid, * FROM jobs WHERE match_score IS NOT NULL ORDER BY rowid DESC")
    results = cursor.fetchall()
    job_id = cursor.lastrowid
    return render_template("jobs.html",results=results, job_id=job_id)


# Retrieve data for specific job ID from list of jobs
@app.route("/jobs/<int:job_id>")
def get_job_by_id(job_id):
    db = get_db()
    cursor = db.cursor() # type: ignore
    cursor.execute("SELECT * FROM jobs WHERE rowid = ?", (job_id,))
    job = cursor.fetchone()

    if job is None:
        flash("Job not found.", "error")
        return redirect(url_for("index"))

    return render_template("submit_form.html", 
        nw_alignment=job["nw_result"],
        nw_seq1=job["nw_seq1_name"],
        nw_seq2=job["nw_seq2_name"],
        nw_score=job["score"],
        nw_percent_id=job["percent_id"],
        nw_gaps=job["gaps"],
        sw_alignment=job["sw_result"],
        sw_score=job["score2"],
        sw_percent_id=job["percent_id2"],
        sw_gaps=job["gaps2"],
        sw_seq1=job["sw_seq1_name"],
        sw_seq2=job["sw_seq2_name"])


# Run algorithm
def main(match, mismatch, gap, sequence1_path, sequence2_path, molecule, is_text=False):
    if molecule == "Protein":
        matrix = get_aa_matrix("resources\\blosum62.txt")
    else:
        matrix = None
    a = AlignNW(molecule, aa_matrix=matrix, match=match, mismatch=mismatch, gap=gap)
    b = AlignSW(molecule, aa_matrix=matrix, match=match, mismatch=mismatch, gap=gap)

    if is_text:
        return "Input Sequence 1", "Input Sequence 2", a.get_alignment(sequence1_path, sequence2_path), \
                "Input Sequence 1", "Input Sequence 2", b.get_alignment(sequence1_path, sequence2_path)    

    sequence1_name, sequence1, is_multiseq = get_fasta_seq(sequence1_path) # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
    sequence2_name, sequence2, is_multiseq = get_fasta_seq(sequence2_path) # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta

    nw = a.get_alignment(sequence1, sequence2)
    sw = b.get_alignment(sequence1, sequence2)
    nw_result = f"{nw[0]}\n{nw[1]}\n{nw[2]}"
    sw_result = f"{sw[0]}\n{sw[1]}\n{sw[2]}"
    score = nw[3]
    percent_id = nw[4]
    gaps = nw[5]
    score2 = sw[3]
    percent_id2 = sw[4]
    gaps2 = sw[5]
    
    return sequence1_name, sequence2_name, nw_result, \
            sequence1_name, sequence2_name, sw_result, \
                score, percent_id, gaps, score2, percent_id2, gaps2


if __name__ == "__main__":
    #main()#match, mismatch, gap)
    app.run(debug=True)