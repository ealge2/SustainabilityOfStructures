# file contains code for generating plot for comparison of cross-sections
# idea of plot: what maximal bending resistance can be reached with a cross-section that uses a budget of 100 kg C02/m2
# output: M-Chi plot with corresponding cross-section geometries

# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import plot_datasets  # file with code for plotting results in a standardized way
import matplotlib.pyplot as plt

#  INPUT
# create dummy-database
database_name = "dummy_sustainability.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

#  define content of plot
gwp_budget = 100  # kg Co2-eq/m2

# DEFINE CONTENT OF PLOT
# define cross-section types
cross_section_types = ["wd_rec", "rc_rec"]

# Create a figure and axis for M-Chi relationship
fig, ax = plt.subplots()

# ADD LABELS, LEGEND AND AXIS LIMITS AND GRID TO THE PLOTS
plt.xlabel('Chi [1/m]')
plt.ylabel("m [kNm/m]")
plt.grid()

# define materials to search for data in the database (table products, attribute material)
for cross_section_type in cross_section_types:
    if cross_section_type[0:2] == "wd":
        mat_names = ["'glue-laminated_timber'", "'solid_structural_timber_(kvh)'"]
        # retrive data from database, find optimal cross-sections and plot results
    elif cross_section_type[0:2] == "rc":
        mat_names = ["'ready_mixed_concrete'"]
    else:
        mat_names = []
        print("relevant materials for cross-section type " + cross_section_type + " are not defined yet")
    plot_datasets.plot_section_dataset(database_name, cross_section_type, mat_names, ax, gwp_budget)

# SHOW FIGURE
plt.show()
