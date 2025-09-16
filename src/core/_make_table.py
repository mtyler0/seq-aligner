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

if __name__ == "__main__":
    create_table()
