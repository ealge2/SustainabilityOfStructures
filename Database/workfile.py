import sqlite3

from ModifyTables import insert_db, product_is_unique

# create or open database sustainability
connection = sqlite3.connect('database_250326.db')
# create cursor object
cursor = connection.cursor()


cursor.execute("DROP TABLE products")

# safe changes in database
connection.commit()

# close database
connection.close()
