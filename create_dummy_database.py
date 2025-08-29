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
    PRO_ID INTEGER PRIMARY KEY, 
    source VARCHAR(30), 
    EPD_date DATE, 
    valid_from DATE, 
    valid_to DATE, 
    product_name VARCHAR(30), 
    oekobaudat_id VARCHAR(6), 
    material_group VARCHAR(20), 
    material VARCHAR(20), 
    cement VARCHAR(20), 
    mech_prop VARCHAR(20), 
    density FLOAT, 
    GWP FLOAT, 
    cost FLOAT,
    cost2 FLOAT);"""
    cursor.execute(sql_command)

    # fill dummy data for concrete C20/25 into table products
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("Ecoinvent", NULL, NULL, NULL, "NULL", "1.4.01", "mineral_building_products",
             "ready_mixed_concrete", "CEM II/A CH-Mix", "C20/25", 2420, 81e-3, 200, 50);"""
    cursor.execute(sql_command)

    # fill dummy data for concrete C25/30 into table products
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("Betonsortenrechner", NULL, NULL, NULL, "NPK B RC-C50", "1.4.01", "mineral_building_products",
             "ready_mixed_concrete", "CEM II/B CH-Mix", "C25/30", 2242, 100e-3, 220, 50);"""
    cursor.execute(sql_command)

    # fill dummy data for concrete C30/37 into table products
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("Ecoinvent", NULL, NULL, NULL, "NULL", "1.4.01", "mineral_building_products",
             "ready_mixed_concrete", "CEM II/B CH-Mix", "C30/37", 2400, 85e-3, 240, 50);"""
    cursor.execute(sql_command)

    # fill dummy data for reinforcing steel B500B into table products
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("KBOB", NULL, NULL, NULL, "Betonstahl KBOB", "4.1.01", "metals", "steel_reinforcing_bar",
             NULL, "B500B", 7850, 773e-3, 11775, NULL);"""
    cursor.execute(sql_command)

    # fill dummy data for timber GL24h into table products (origin CH)
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("KBOB", NULL, NULL, NULL, "Brettschichtholz KBOB, CH", "3.1.04", "wood",
            "glue-laminated_timber", NULL, "GL24h", 439, 253e-3, 1200, 15);"""
    cursor.execute(sql_command)

    # fill dummy data for timber GL24h into table products
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("KBOB", NULL, NULL, NULL, "Brettschichtholz KBOB", "3.1.04", "wood", "glue-laminated_timber",
             NULL, "GL24h", 439, 335e-3, 1100, 15);"""
    cursor.execute(sql_command)

    # fill dummy data for timber C24 into table products
    sql_command = """INSERT INTO products (source, EPD_date, valid_from, valid_to, product_name, oekobaudat_id,
            material_group, material, cement, mech_prop, density, GWP, cost, cost2)
            VALUES ("KBOB", NULL, NULL, NULL, "Konstruktionsvollholz KBOB", "3.1.02", "wood",
             "solid_structural_timber_(kvh)", NULL, "C24", 436, 156e-3, 1000, 15);"""
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

    # fill data for concrete C20/25 into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("C20/25", 20e6 , 2.3e6, NULL,
        NULL, 30e9, 25e3, NULL);"""
    cursor.execute(sql_command)

    # fill data for concrete C25/30 into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("C25/30", 25e6 , 2.6e6, NULL,
        NULL, 32e9, 25e3, NULL);"""
    cursor.execute(sql_command)

    # fill data for concrete C30/37 into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("C30/37", 25e6 , 2.9e6, NULL,
        NULL, 34e9, 25e3, NULL);"""
    cursor.execute(sql_command)

    # fill data for B500B into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("B500B", 500e6 , 500e6, NULL,
        NULL, 205e9, NULL, NULL);"""
    cursor.execute(sql_command)

    # fill data for timber C16 into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("C16", NULL , NULL, 16e6,
        1.5e6, 8e9, 5e3, 0.8e-3);"""
    cursor.execute(sql_command)

    # fill data for timber C24 into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("C24", NULL , NULL, 24e6,
        1.5e6, 11e9, 5e3, 0.8e-3);"""
    cursor.execute(sql_command)

    # fill data for timber GL24h into table material_prop
    sql_command = """INSERT INTO material_prop (name, strength_comp, strength_tens, strength_bend,
     strength_shea, E_modulus, density_load, burn_rate)
        VALUES ("GL24h", NULL , NULL, 24e6,
        1.8e6, 11.5e9, 5e3, 0.7e-3);"""
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
