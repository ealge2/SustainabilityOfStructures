# file contains code for generating fundamental plots of first example (simple supported beam,
# reinforced concrete cross-section)

# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import plot_datasets  # file with code for plotting results in a standardized way
import matplotlib.pyplot as plt

# INPUT
# create dummy-database
database_name = "dummy_sustainability_1.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# define system lengths for plot
lengths = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

# max. number of iterations per optimization. Higher value leads to better results
max_iter = 20

#  define content of plot
criteria = ["ENV"]
optima = ["GWP"]

# create floor structure for solid wooden cross-section
bodenaufbau_rc_slap = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                       ["'Unterlagsboden Zement, 85 mm'", False, False], ["'Glaswolle'", 0.03, False]]
bodenaufbau_rc = struct_analysis.FloorStruc(bodenaufbau_rc_slap, database_name)

# define loads on member
g2k = 0.75e3  # n.t. Einbauten
qk = 2e3  # Nutzlast

# define service limit state criteria
req = struct_analysis.Requirements()

# define materials for which date is searched in the database (table products, attribute material)
mat_names = ["'ready_mixed_concrete'"]
# retrieve data from database, find optimal cross-sections and plot results
data_max, vrfctn_members = plot_datasets.plot_dataset(lengths, database_name, criteria, optima, bodenaufbau_rc, req,
                                                      "rc_rec", mat_names, g2k, qk, max_iter)

# define legend of plots
plotted_data = [["h_struct", "[m]"], ["h_tot", "[m]"], ["GWP_struct", "[kg-CO2-eq]"], ["GWP_tot", "[kg-CO2-eq]"],
                ["cost_struct", "[CHF]"]]

for idx, info in enumerate(plotted_data):
    plt.subplot(3, 2, idx + 1)
    plt.xlabel('l [m]')
#    plt.title(info[0])
    plt.ylabel(info[0] + " " + info[1])
    if idx % 2 == 0:
        plt.axis((min(lengths), max(lengths), 0, max(data_max[idx], data_max[idx+1])))
    else:
        plt.axis((min(lengths), max(lengths), 0, max(data_max[idx], data_max[idx-1])))
    # plt.legend()
    plt.grid()

# plot cross-section of members for verification
for i, mem in enumerate(vrfctn_members[0]):
    section = mem.section
    plot_datasets.plot_section(section)
    # Show the plot
    plt.title(f'#{vrfctn_members[1][i]}')

plt.show()
