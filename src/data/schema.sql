CREATE TABLE alignments (
    rowid INTEGER PRIMARY KEY,
    current_time TEXT DEFAULT (datetime('now', '-5 hours')),
    nw_seq1_name TEXT,
    nw_seq2_name TEXT,
    nw_result TEXT,
    sw_seq1_name TEXT,
    sw_seq2_name TEXT,
    sw_result TEXT
);