
from create_dummy_database import create_database

create_database('test.db')


"""import sqlite3

connection = sqlite3.connect("test.db")
cursor = connection.cursor()

for row in cursor.execute('SELECT * FROM products;'):
    print(row)"""


