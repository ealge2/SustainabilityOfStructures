# file creates a dummy database for testing the structure analysis code
# units: [m], [kg], [s], [N], [CHF]
import sqlite3


def create_database_slab(database_name):
    # create or open database sustainability
    connection = sqlite3.connect(database_name)

    # create cursor object
    cursor = connection.cursor()

    # delete existing products table
    try:
        cursor.execute("""DROP TABLE slab_properties;""")
    except:
        pass

    # create table for products data
    sql_command = """
    CREATE TABLE slab_properties ( 
    NAME TEXT, 
    RAENDER TEXT, 
    LX FLOAT, 
    LY FLOAT, 
    MX_POS FLOAT,
    MY_POS FLOAT,
    MX_NEG FLOAT,
    MY_NEG FLOAT,
    V_POS FLOAT,
    V_NEG FLOAT,
    W FLOAT,
    F FLOAT);"""
    cursor.execute(sql_command)

    # fill slab properties into db

    entries = [
        ("4S_3x3", "LL-frei", 3, 3, 0.04244, 0, 0.04244, 0, 0.30633, 0.30633, 0.00394, 0),
        ("4S_5x5", "LL-frei", 5, 5, 0.04296, 0, 0.04296, 0, 0.3206, 0.30633, 0.00391, 0),
        ("4S_6x6", "LL-frei", 6, 6, 0.042944, 0, 0.042944, 0, 0.3235, 0.3235, 0.00345, 0),
        ("4S_7x7", "LL-frei", 7, 7, 0.04298, 0, 0.04298, 0, 0.326, 0.326, 0.003945, 0),
        ("4S_8x8", "LL-frei", 8, 8, 0.043047, 0, 0.043047, 0, 0.3275, 0.3275, 0.00393, 0),
        ("4S_10x10", "LL-frei", 10, 10, 0.04304, 0, 0.04304, 0, 0.3297, 0.3297, 0.003945, 0)
    ]

    for entry in entries:
        sql_command = """INSERT INTO slab_properties (NAME, RAENDER, LX, LY, MX_POS, MY_POS, MX_NEG, MY_NEG, V_POS, V_NEG, W, F )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"""
        cursor.execute(sql_command, entry)

    # safe changes in database
    connection.commit()

    # close database
    connection.close()


database_name = "slab_properties.db"
create_database_slab(database_name)



import sqlite3
from tabulate import tabulate  # pip install tabulate (optional, für schöne Ausgabe)

def show_database_contents(database_name):
    connection = sqlite3.connect(database_name)
    cursor = connection.cursor()


    cursor.execute("SELECT * FROM slab_properties")
    rows = cursor.fetchall()

    # Spaltennamen holen
    column_names = [description[0] for description in cursor.description]

    # Ausgabe als Tabelle
    print(tabulate(rows, headers=column_names, tablefmt="grid"))

    connection.close()

# Aufruf
show_database_contents("slab_properties.db")
