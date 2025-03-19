import sqlite3

# create or open database sustainability
connection = sqlite3.connect('data.db')
# create cursor object
cursor = connection.cursor()


#checks if Product ID is unique
def product_is_unique(index):
    rows = cursor.execute("SELECT * FROM products").fetchall()
    print("Total rows are:  ", len(rows))
    for product in rows:
        print("2")
        if product[0] == index:
            return False
        return True

def insert_db():
    index = input("Index >>")

    if product_is_unique(index):
        EPD_ID = input("EPD_ID >")
        source = input("Source >")
    else:
        print("product already exists")


 #   connection.commit()

def edit_db():
        EPD_ID = input("Which field would you like to edit?")

def get_product_info_db():
    target_product = input("What product do you want to see information about? >>")

def display_db():
    rows = cursor.execute("SELECT * FROM products").fetchall()
    print("Products: ")
    for user in rows:
        print(f"-{user[0]} - {user[1]}- {user[2]}")


def exit_db():
    cursor.close()
    connection.close()


def select_options():
        options = input("""
        -------------------------------------
        Type '0' to exit
        Type '1' to insert a new product
        Type '2' to display products
        Type '3' to delete product
        Type '4' to edit product
        Type '5' to get product information
        -------------------------------------
        >>""")