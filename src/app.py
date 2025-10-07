from core.utils import *
from core.needleman_wunsch import *
from core.smith_waterman import *
from flask import Flask, render_template, request, redirect, g, flash, url_for
from psycopg.rows import dict_row


app = Flask(__name__)
app.secret_key = "2nf38-f2n398"
app.config["UPLOAD_FOLDER"] = "data/uploaded_files"


@app.route("/", methods=["GET", "POST"])
def index():
    return render_template("index.html")


# Send data to algorithm, store in db, then redirect to submit page
@app.route("/submit_form", methods=["POST"])
def post_form():
    upload_folder = app.config["UPLOAD_FOLDER"]
    text_input: str = request.form["SUBJECT"].strip() and request.form["QUERY"].strip()
    input: bool = len(text_input) > 0 or (
        request.files.get("seq1file").filename != "" and request.files.get("seq2file").filename != "" #type: ignore
        ) 

    # Scoring params
    match = request.form.get("match", type=int)
    mismatch = request.form.get("mismatch", type=int)
    gap = request.form.get("gap", type=int)
    molecule = request.form.get("type")
    matrix = request.form.get("protein-matrix")

    if not input:
        flash("Error: Missing file/text input")
        return redirect("/")

    try:
        if text_input:
            seq1_text = request.form["SUBJECT"]
            seq2_text = request.form["QUERY"]
            params = main(match, mismatch, gap, seq1_text, seq2_text, molecule, is_text=True)
        else:
            seq1_file = request.files["seq1file"]
            seq2_file = request.files["seq2file"]
            seq1_path = save_file(seq1_file, upload_folder)
            seq2_path = save_file(seq2_file, upload_folder)
            params = main(match, mismatch, gap, seq1_path, seq2_path, molecule)

    # Use nested with statements here
        with get_db() as conn:
            with conn.cursor(row_factory=dict_row) as c:
                c.execute(
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
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) RETURNING id;
                    """, (match, mismatch, gap, *params)
                    )
                job_id = c.fetchone()["id"] # type: ignore

        return redirect(url_for("submit", job_id=job_id))

    except Exception as e:
        print(f"Type of e: {type(e)}, value: {e}")
        flash(f"ERROR: {e}")
        return redirect(url_for("index"))


# Get alignment data for the current job ID
@app.route("/submit")
def submit():
    job_id = request.args.get("job_id")
    print(job_id)
    try:
        with get_db() as conn:
            with conn.cursor(row_factory=dict_row) as c:
                c.execute("SELECT * FROM jobs WHERE id = %s", (job_id,))
                job = c.fetchone()

        if job is None:
            flash("Job not found.", "error")
            return redirect(url_for("index"))

        return render_template(
            "submit_form.html", 
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
            sw_seq2=job["sw_seq2_name"]
            )
        
    except Exception as e:
        print(f"ERROR: {e}")
        flash(f"ERROR: {e}")
        return redirect(url_for("index"))


# Retrieve all past data
@app.route("/jobs")
def get_jobs():
    try:
        with get_db() as conn:
            with conn.cursor(row_factory=dict_row) as c:
                c.execute("SELECT id, * FROM jobs WHERE match_score IS NOT NULL ORDER BY id DESC")
                results = c.fetchall()

        return render_template("jobs.html",results=results)

    except Exception as e:
        print(f"ERROR: {e}")
        flash(f"ERROR: {e}")
        return redirect(url_for("index"))


# Retrieve data for specific job ID from list of jobs
@app.route("/jobs/<int:job_id>")
def get_job_by_id(job_id):
    try:
        with get_db() as conn:
            with conn.cursor(row_factory=dict_row) as c:
                c.execute("SELECT * FROM jobs WHERE id = %s", (job_id,))
                job = c.fetchone() # type: ignore

        if job is None:
            flash("Job not found.", "error")
            return redirect(url_for("index"))

        return render_template(
            "submit_form.html", 
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
            sw_seq2=job["sw_seq2_name"]
            )
        
    except Exception as e:
        print(f"ERROR: {e}")
        flash(f"ERROR: {e}")
        return redirect(url_for("index"))


@app.teardown_appcontext
def close_conneciton(exception):
    db = getattr(g, "_database", None)
    if db is not None:
        db.close()


# Run algorithm
def main(match, mismatch, gap, sequence1_path, sequence2_path, molecule, is_text=False):
    if molecule == "Protein":
        matrix = get_aa_matrix("resources\\blosum62.txt")
    else:
        matrix = None

    a = AlignNW(molecule, aa_matrix=matrix, match=match, mismatch=mismatch, gap=gap)
    b = AlignSW(molecule, aa_matrix=matrix, match=match, mismatch=mismatch, gap=gap)

    if is_text:
        sequence1_name = "Input Sequence 1"
        sequence2_name = "Input Sequence 2"
        sequence1 = sequence1_path
        sequence2 = sequence2_path
    else:
        sequence1_name, sequence1 = get_fasta_seq(sequence1_path) # DNA1: data\\seq1.fasta, Protein: data\\human_hbb.fasta
        sequence2_name, sequence2 = get_fasta_seq(sequence2_path) # DNA2: data\\seq2.fasta, Protein: data\\puffer_hbb.fasta

    nw = a.get_alignment(sequence1, sequence2)
    sw = b.get_alignment(sequence1, sequence2)
    nw_result = nw[0]
    sw_result = sw[0]
    score = nw[1]
    percent_id = nw[2]
    gaps = nw[3]
    score2 = sw[1]
    percent_id2 = sw[2]
    gaps2 = sw[3]
    
    return sequence1_name, sequence2_name, nw_result, \
            sequence1_name, sequence2_name, sw_result, \
                score, percent_id, gaps, score2, percent_id2, gaps2


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)