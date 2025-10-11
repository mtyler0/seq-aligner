# Sequence Alignment Web Application

A lightweight, full-stack Flask web app that performs **pairwise sequence alignments** between nucleotide or protein sequences using a dynamic programming algorithm. The application allows users to input FASTA sequences, select alignment parameters, and view formatted alignment results directly in the browser. Designed with modularity, clean backend logic, and scalability in mind.

---

## 🚀 Features

- **Interactive Alignment Form:** Upload or paste sequences, select target database and scoring matrix dynamically.
- **Backend Integration:** Uses implementations of the Needleman-Wunsch (global) and Smith-Waterman (local) algorithms for performing optimal sequence alignments.
- **Customizable Parameters:** Modify match, mismatch, and gap penalties and choose between `BLOSUM62` and `PAM160` substitution matrices for protein scoring
- **Robust Error Handling:** Gracefully manages invalid inputs, failed alignments, and improper file formats.
- **Session-Based Job Tracking:** Users can view and reload past alignments during a session.
- **Responsive Frontend:** Built with templated HTML and CSS for a seamless user experience.
- **Maintainable Architecture:** Modular Flask structure with separated concerns (`index`, `submit`, `jobs`, `get_job_by_id`).
- **Containerized Deployment:** Fully Dockerized for reproducibility and portability.

---

## 🧠 Tech Stack

| Layer | Technology |
|-------|-------------|
| **Frontend** | HTML5, CSS3, Jinja2 |
| **Backend** | Python 3.13, Flask |
| **Database** | PostgreSQL (Session Tracking) |
| **Testing** | Pytest |
| **Version Control** | Git, GitHub |

---

## 📁 Project Structure

```
seq-aligner/
├── init/
│ └── schema.sql
│
├── src/
│ ├── app.py
│ │
│ ├── core/
│ │ ├── migrate.py
│ │ ├── aligner_base.py
│ │ ├── needleman_wunsch.py
│ │ ├── smith_waterman.py
│ │ └── utils.py
│ │
│ ├── data/
│ │ ├── examples/
│ │ ├── uploaded_files/
│ │ └── saved_jobs.db
│ │
│ ├── resources/
│ │ ├── blosum62.txt
│ │ └── pam160.txt
│ │
│ ├── static/
│ │ ├── images/
│ │ └── styles.css
│ │
│ └── templates/
│ ├── base.html
│ ├── index.html
│ ├── jobs.html
│ └── submit_form.html
│
├── tests/
│ └── test_align.py
│
├── .gitignore
├── docker-compose.yml
├── Dockerfile
├── README.md
└── requirements.txt
```

---

## ⚙️ Installation

1. **Clone the repository**
```bash
git clone https://github.com/<your-username>/sequence-aligner.git
cd sequence-aligner
```
2. **Build and run with docker**
```bash
docker compose --build -d
```
3. Access the app at http://localhost:5000

---

## 🧩 Example Workflow

1. **User submits input sequences and parameters**
   - Enter or upload query and subject sequences.
2. **Backend validates and processes request**
   - Input is sanitized and passed to the `get_alignment()` function in `/core`.
3. **Alignment performed using the AlignNW/AlignSW classes**
   - Results are formatted and returned as HTML or error messages.
4. **Results displayed**
   - Users can view formatted alignments and revisit past results during session.

---

## 🔍 Error Handling

- `get_alignment()` raises `ValueError` or `KeyError` for invalid or failed alignments.
- `main()` catches and re-raises descriptive exceptions.
- Flask `submit_form` route flashes user-friendly error messages.
- Ensures clean separation of concerns and readable tracebacks for debugging.

---

## 🧪 Testing

Unit tests are implemented using **pytest** to validate:
- Input validation
- Accurate alignment scoring
- Error propagation and flash messaging

Run tests:
```bash
pytest tests\test_align.py
```

---

## 🔧 Future Improvements

- User authentication and session management
- REST API for remote alignment submission
- Enhance alignment visualization (e.g., color-coded residue matches)
- Integration with cloud storage for file handling
- Implement caching layer for repeated queries
- Deployment on AWS / Render / DigitalOcean

---

## 👤 Author

**Michael Tyler**  
Biochemist turned developer building tools that bridge molecular science and software engineering.  
🔗 [LinkedIn](https://linkedin.com/in/mtyler0)  
💻 [GitHub](https://github.com/mtyler0)

---

## 🧾 License

This project is licensed under the **MIT License**.
