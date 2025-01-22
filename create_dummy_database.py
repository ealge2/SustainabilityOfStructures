# file creates a dummy database for testing the structure analysis code
# units: [m], [kg], [s], [N], [CHF]
import sqlite3

def create_database(data_base_name):
    # create or open database sustainability
    connection = sqlite3.connect(data_base_name)

    # create cursor object
    cursor = connection.cursor()

    # delete existing products table
    try:
        cursor.execute("""DROP TABLE products;""")
    except:
        pass

    # create table for products data
    sql_command = """
    CREATE TABLE products ( 
    EPD_ID INTEGER PRIMARY KEY, 
    source VARCHAR(30), 
    EPD_date DATE, 
    valid_from DATE, 
    valid_to DATE, 
    product_name VARCHAR(30), 
    material VARCHAR(20), 
    kind VARCHAR(20), 
    cement VARCHAR(20), 
    mech_prop VARCHAR(20), 
    density FLOAT, 
    GWP FLOAT, 
    cost FLOAT,
    cost2 FLOAT);"""
    cursor.execute(sql_command)

    # fill dummy data for concrete C25/30 into table products
    sql_command = """INSERT INTO products (EPD_ID, source, EPD_date, valid_from, valid_to, product_name, material, kind,
         cement, mech_prop, density, GWP, cost, cost2)
            VALUES (NULL, "Betonsortenrechner", NULL, NULL, NULL, "NPK B RC-C50", "concrete", "structural",
             "CEM II/B CH-Mix", "C25/30", 2190, 98e-3, 220, 50);"""
    cursor.execute(sql_command)

    # fill dummy data for reinforcing steel B500B into table products
    sql_command = """INSERT INTO products (EPD_ID, source, EPD_date, valid_from, valid_to, product_name, material, kind,
         cement, mech_prop, density, GWP, cost, cost2)
            VALUES (NULL, "KBOB", NULL, NULL, NULL, "Betonstahl KBOB", "metal", "reinforcing steel",
             NULL, "B500B", 7850, 773e-3, 11775, NULL);"""
    cursor.execute(sql_command)

    # fill dummy data for timber GL24h into table products
    sql_command = """INSERT INTO products (EPD_ID, source, EPD_date, valid_from, valid_to, product_name, material, kind,
         cement, mech_prop, density, GWP, cost, cost2)
            VALUES (NULL, "KBOB", NULL, NULL, NULL, "Brettschichtholz KBOB", "timber", "glulam",
             NULL, "GL24h", 439, 253e-3, 1200, 15);"""
    cursor.execute(sql_command)

    try:
        # delete existing material properties table
        cursor.execute("""DROP TABLE material_prop;""")
    except:
        pass

    # create table for material properties data
    sql_command = """
    CREATE TABLE material_prop ( 
    Mat_ID INTEGER PRIMARY KEY, 
    name VARCHAR(20), 
    strength_comp FLOAT, 
    strength_tens FLOAT,
    strength_bend FLOAT, 
    strength_shea FLOAT,
    E_modulus FLOAT,
    density_load FLOAT,
    burn_rate FLOAT);"""
    cursor.execute(sql_command)

    # fill data for concrete C25/30 into table material_prop
    sql_command = """INSERT INTO material_prop (Mat_ID, name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES (NULL, "C25/30", 25e6 , 2.6e6, NULL,
        NULL, 30e9, 25e3, NULL);"""
    cursor.execute(sql_command)

    # fill data for B500B into table material_prop
    sql_command = """INSERT INTO material_prop (Mat_ID, name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES (NULL, "B500B", 500e6 , 500e6, NULL,
        NULL, 205e9, NULL, NULL);"""
    cursor.execute(sql_command)

    # fill data for timber GL24h into table material_prop
    sql_command = """INSERT INTO material_prop (Mat_ID, name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES (NULL, "GL24h", NULL , NULL, 24e6,
        1.8e6, 11e9, 5e3, 0.7e-3);"""
    cursor.execute(sql_command)

    # delete existing floor structure property table
    try:
        cursor.execute("""DROP TABLE floor_struc_prop;""")
    except:
        pass

    # create table for floor structure materials data
    sql_command = """
    CREATE TABLE floor_struc_prop ( 
    Mat_ID INTEGER PRIMARY KEY, 
    name VARCHAR(20), 
    h_fix FLOAT,
    E FLOAT,
    density FLOAT,
    weight,
    GWP FLOAT);"""
    cursor.execute(sql_command)

    # fill data for parquet into table floor_struc_mat
    sql_command = """INSERT INTO floor_struc_prop (Mat_ID, name, h_fix, E, density, weight, GWP)
        VALUES (NULL, "Parkett 2-Schicht werkversiegelt, 11 mm", 0.011, NULL, 555 , 8e3, 1279e-3);"""
    cursor.execute(sql_command)

    # fill data for screed into table floor_struc_mat
    sql_command = """INSERT INTO floor_struc_prop (Mat_ID, name, h_fix, E, density, weight, GWP)
        VALUES (NULL, "Unterlagsboden Zement, 85 mm", 0.085, 21e9, 1850, 22e3, 120e-3);"""
    cursor.execute(sql_command)

    # fill data for impact sound insulation into table floor_struc_mat
    sql_command = """INSERT INTO floor_struc_prop (Mat_ID, name, h_fix, E, density, weight, GWP)
        VALUES (NULL, "Glaswolle", NULL, NULL, 80, 0.8e3, 1100e-3);"""
    cursor.execute(sql_command)

    # fill data for grit into table floor_struc_mat
    sql_command = """INSERT INTO floor_struc_prop (Mat_ID, name, h_fix, E, density, weight, GWP)
        VALUES (NULL, "Kies gebrochen", NULL, NULL, 2000, 20e3, 18e-3);"""
    cursor.execute(sql_command)

    # safe changes in database
    connection.commit()

    # close database
    connection.close()
