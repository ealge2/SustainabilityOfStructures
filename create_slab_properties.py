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

    # fill dummy data for concrete C20/25 into table products
    sql_command = """INSERT INTO slab_properties (NAME, RAENDER, LX, LY, MX_POS, MY_POS, MX_NEG, MY_NEG, V_POS, V_NEG, W, F )
            VALUES ("4S_3x3", "LL-frei", 3, 3, 0.04244, 0, 0.04244, 0, 0.30633, 0.30633, 0.00394, 0);"""
    cursor.execute(sql_command)



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
