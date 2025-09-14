CREATE TABLE jobs (
    current_time TEXT DEFAULT (datetime('now', 'localtime')),
    match_score INT,
    mismatch INT,
    gap INT,
    nw_seq1_name TEXT,
    nw_seq2_name TEXT,
    nw_result TEXT,
    sw_seq1_name TEXT,
    sw_seq2_name TEXT,
    sw_result TEXT
);
