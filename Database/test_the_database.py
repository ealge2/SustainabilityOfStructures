import sqlite3  # import modul for SQLite

database_name = "database_250326.db"

connection = sqlite3.connect(database_name)
cursor = connection.cursor()


# retrieve column names of table products
table_name = "products"
cursor.execute(f"PRAGMA table_info({table_name})")
columns = [row[1] for row in cursor.fetchall()]
print("Columns:", columns)

## ToDo: change culumn names of database to the ones defined in excel. the format is not part of the column name!


# test how to get values for a specific steel
mat_name = "'B500B'"
inquiry = ("SELECT country FROM products WHERE"
           " mech_prop=" + mat_name)
cursor.execute(inquiry)
result = cursor.fetchall()
print(result)
# -> works

# test how to get values for a specific material
mat_name = "'B500B'"
inquiry = ("SELECT country FROM products WHERE"
           " mech_prop=" + mat_name)
cursor.execute(inquiry)
result = cursor.fetchall()
print(result)
# -> works



