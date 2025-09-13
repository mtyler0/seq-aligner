import sqlite3

def create_table():
    # Connect (creates the DB file if it doesn’t exist)
    conn = sqlite3.connect("data\\saved_jobs.db")
    cursor = conn.cursor()

    # Create the table if it doesn’t exist
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            current_time TEXT DEFAULT (datetime('now', 'localtime')),
            nw_seq1_name TEXT,
            nw_seq2_name TEXT,
            nw_result TEXT,
            sw_seq1_name TEXT,
            sw_seq2_name TEXT,
            sw_result TEXT
        )
    """)

    conn.commit()
    conn.close()
    print("Table 'jobs' created (if it didn't already exist).")

if __name__ == "__main__":
    create_table()
