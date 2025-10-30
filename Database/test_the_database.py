import sqlite3  # import modul for SQLite

database_name = 'database_251030.db'

connection = sqlite3.connect(database_name)
cursor = connection.cursor()

# retrieve column names of table products
table_name = "products"
cursor.execute(f"PRAGMA table_info({table_name})")
columns = [row[1] for row in cursor.fetchall()]
print("Columns:", columns)

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

# Referenz aus Code Luisa (funktioniert)
KBOB_steel = cursor.execute("""
                        SELECT A1toA3_GWP FROM products
                        WHERE MATERIAL LIKE '%structural steel profile%'
                        AND "source [string]" LIKE '%KBOB%'
                        """).fetchall()
print(KBOB_steel)

#------------------------------------------------------------------------------------------------------------
# Abwandlung für Integration in eigenen Code (funktioniert)
inquiry = ("""
        SELECT A1toA3_GWP FROM products
        WHERE  "material [string]" LIKE 'ready_mixed_concrete'
        """)
KBOB_steel = cursor.execute(inquiry).fetchall()
print(KBOB_steel)

# Abwandlung für Integration in eigenen Code (zu testen)
mat_name = "'ready_mixed_concrete'"
inquiry = ("""
        SELECT A1toA3_GWP FROM products
        WHERE  "material [string]" LIKE """ + mat_name
        )
KBOB_steel = cursor.execute(inquiry).fetchall()
print(KBOB_steel)

# Abwandlung für Integration in eigenen Code (zu testen)
mat_name = "'glue_laminated_timber'"
inquiry = ("""
        SELECT PRO_ID FROM products
        WHERE  "material [string]" LIKE """ + mat_name
        )
KBOB_steel = cursor.execute(inquiry).fetchall()
print(KBOB_steel)


