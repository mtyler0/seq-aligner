import sqlite3

def create_table():
    # Connect (creates the DB file if it doesn’t exist)
    conn = sqlite3.connect("C:\\Users\\mtyle\\Desktop\\ccode\\NW_algo\\src\\data\\saved_jobs.db")
    cursor = conn.cursor()

    # Create the table if it doesn’t exist
    with open("C:\\Users\\mtyle\\Desktop\\ccode\\NW_algo\\src\\data\\schema.sql", "r") as file:
        r = file.read()
        cursor.executescript(r)
        conn.commit()
        conn.close()

    print("Table 'jobs' created (if it didn't already exist).")


def add_column():
    conn = sqlite3.connect("C:\\Users\\mtyle\\Desktop\\ccode\\NW_algo\\src\\data\\saved_jobs.db")
    cursor = conn.cursor()
    
    columns = [("score2", "INT"), ("percent_id2", "TEXT"), ("gaps2", "INT")]

    for column, dtype in columns:
        try:
            cursor.execute(f"ALTER TABLE jobs ADD COLUMN {column} {dtype};")
            print(f"Column '{column}' added successfully.")
        except sqlite3.OperationalError as e:
            print(f"Failed on {column}: {e}")
            
    conn.commit()
    conn.close()

    
if __name__ == "__main__":
    pass
