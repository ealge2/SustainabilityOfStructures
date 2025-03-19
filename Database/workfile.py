import sqlite3

from ModifyTables import insert_db, product_is_unique

# create or open database sustainability
connection = sqlite3.connect('data.db')
# create cursor object
cursor = connection.cursor()


rows = cursor.execute("SELECT * FROM products").fetchall()

print("Total rows are:  ", len(rows))

for name in rows:
    print("index: ", name[0])
    if name[0] == 52:
        print("false")


product_is_unique(2)

