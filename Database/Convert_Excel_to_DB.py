#Packages: pip install pandas openpyxl

import pandas as pd
import sqlite3


# Pfad zur Excel-Datei (bitte anpassen!)
excel_file = "Datenbankdefinition_lua_250423.xlsx"   #name of excel file. must be in same directory as code

# Excel-Datei einlesen
df_products = pd.read_excel(excel_file, sheet_name="products", engine="openpyxl")
df_material_prop = pd.read_excel(excel_file, sheet_name="material_prop", engine="openpyxl")

# Verbindung zur SQLite-Datenbank herstellen (oder neu erstellen)
conn = sqlite3.connect("database.db") #name of database

# Daten in die Datenbank schreiben
df_products.to_sql("products", conn, if_exists="replace", index=False)
df_material_prop.to_sql("material_prop", conn, if_exists="replace", index=False)

# Verbindung schließen
conn.close()

print("Datenbank erfolgreich erstellt")

