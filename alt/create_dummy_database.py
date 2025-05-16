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
    Cost FLOAT);"""
    cursor.execute(sql_command)

    # fill dummy data for concrete C30/37 into table products
    sql_command = """INSERT INTO products (EPD_ID, source, EPD_date, valid_from, valid_to, product_name, material, kind,
     cement, mech_prop, density, GWP, Cost)
        VALUES (NULL, "Betonsortenrechner", NULL, NULL, NULL, "NPK B RC-C50", "concrete", "structural",
         "CEM II/B CH-Mix", "C25/30", 2190, 98, 220);"""
    cursor.execute(sql_command)

    # fill dummy data for reinforcing steel B500B into table products
    sql_command = """INSERT INTO products (EPD_ID, source, EPD_date, valid_from, valid_to, product_name, material, kind,
     cement, mech_prop, density, GWP, Cost)
        VALUES (NULL, "KBOB", NULL, NULL, NULL, "Betonstahl KBOB", "metal", "reinforcing steel",
         NULL, "B500B", 7850, 773, 11775);"""
    cursor.execute(sql_command)

    # fill dummy data for timber GL24h into table products
    sql_command = """INSERT INTO products (EPD_ID, source, EPD_date, valid_from, valid_to, product_name, material, kind,
     cement, mech_prop, density, GWP, Cost)
        VALUES (NULL, "KBOB", NULL, NULL, NULL, "Brettschichtholz KBOB", "timber", "glulam",
         NULL, "GL24h", 439, 253, 1000);"""
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
    density_load FLOAT);"""
    cursor.execute(sql_command)

    # fill data for concrete C30/37 into table material_prop
    sql_command = """INSERT INTO material_prop (Mat_ID, name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load)
        VALUES (NULL, "C25/30", 25 , 2.6, NULL,
        NULL, 30000, 25);"""
    cursor.execute(sql_command)

    # fill data for B500B into table material_prop
    sql_command = """INSERT INTO material_prop (Mat_ID, name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load)
        VALUES (NULL, "B500B", 500 , 500, NULL,
        NULL, 205000, NULL);"""
    cursor.execute(sql_command)

    # fill data for timber GL24h into table material_prop
    sql_command = """INSERT INTO material_prop (Mat_ID, name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load)
        VALUES (NULL, "GL24h", NULL , NULL, 24,
        1.8, 11000, 5);"""
    cursor.execute(sql_command)



    # safe changes in database
    connection.commit()

    # close database
    connection.close()
