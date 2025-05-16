#Data analysis tool
import pandas as pd
#Library that provieds disk-based database that does not require a separate server process
import sqlite3

#Takes data from csv and moves it into a Pandas dataframe
print("Converting csv to Dataframe...")
df = pd.read_csv(
    filepath_or_buffer='products.csv',
    sep= ';',
    encoding= 'latin1',
    header=0
)

#Inserts the Pandas Dataframe values into SQLITE
print("Connecting to sqlite datatbase...")
connection = sqlite3.connect('data.db')

df.to_sql(
    name = 'products',
    con= connection,
    if_exists = 'replace' #'replace' replaces entire table,'append' checks for duplicate rows and removes them
)

# safe changes in database
connection.commit()

# close database
connection.close()