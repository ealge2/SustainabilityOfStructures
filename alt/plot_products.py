import matplotlib.pyplot as plt
import sqlite3
import numpy as np
import statistics

# create or open database sustainability
#connection = sqlite3.connect('data.db')
connection = sqlite3.connect('../Database/database_250326.db')
# create cursor object
cursor = connection.cursor()

#------------------------------------------------------------------------------------------------------------------------
#extract values for concrete

emissions_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material_0 LIKE '%concrete%' "
                    "AND A1toA3_GWP IS NOT NULL "
                    "AND A1toA3_GWP != 0 "
                    "AND type LIKE '%ready mixed concrete%' "
                    ).fetchall()
emissions_concrete_values = [row[0] for row in emissions_concrete]

EPD_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material_0 LIKE '%concrete%' "
                    "AND A1toA3_GWP IS NOT NULL "
                    "AND A1toA3_GWP != 0 "
                    "AND 'source [string]' NOT LIKE '%Betonsortenrechner%' "
                    "AND 'source [string]' NOT LIKE '%Ecoinvent%' "
                    "AND 'source [string]' NOT LIKE '%KBOB%' "
                    "AND type LIKE '%ready mixed concrete%' "
                    ).fetchall()
EPD_concrete_values = [row[0] for row in EPD_concrete]

KBOB_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material_0 LIKE '%concrete%' "
                    "AND 'source [string]' LIKE '%KBOB%' "
                    ).fetchall()
KBOB_concrete_values = [row[0] for row in KBOB_concrete]

Ecoinvent_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material_0 LIKE '%concrete%' "
                    "AND 'source [string]' LIKE '%Ecoinvent%' "
                    "AND A1toA3_GWP != 0 "
                    ).fetchall()
Ecoinvent_concrete_values = [row[0] for row in Ecoinvent_concrete]

Betonsortenrechner_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material_0 LIKE '%concrete%' "
                    "AND A1toA3_GWP != 0 "
                    "AND 'source [string]' LIKE '%Betonsortenrechner%' "
                    ).fetchall()
Betonsortenrechner_concrete_values = [row[0] for row in Betonsortenrechner_concrete]

#------------------------------------------------------------------------------------------------------------------------
#extract values for wood
emissions_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    "WHERE material_0 LIKE '%timber%' "
                    "AND Total_GWP IS  NOT NULL "
                    "AND 'source [string]' NOT LIKE '%Studiengemeinschaft%' "
                    ).fetchall()
emissions_timber_values = [row[0] for row in emissions_timber]

EPD_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    " WHERE material_0 LIKE '%timber%' "
                    "AND Total_GWP IS NOT NULL "
                    "AND 'source [string]' NOT LIKE '%KBOB%' "
                    "AND EPD_ID NOT LIKE '%verifizierung%' "
                    ).fetchall()
EPD_timber_values = [row[0] for row in EPD_timber]

KBOB_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    "WHERE material_0 LIKE '%timber%' "
                    "AND 'source [string]' LIKE '%KBOB%' "
                    ).fetchall()
KBOB_timber_values = [row[0] for row in KBOB_timber]

print(len(emissions_timber_values))
print(len(EPD_timber_values))
print(emissions_timber_values)
#------------------------------------------------------------------------------------------------------------------------
#extract values for reinforcement

emissions_reinf = cursor.execute("SELECT A1toA3_GWP FROM products "
                           " WHERE 'product_name [string]' LIKE '%Betonstahl%' "
                           "AND A1toA3_GWP IS NOT NULL "
                           "AND A1toA3_GWP != 0 "
                           ).fetchall()
emissions_reinf_values = [row[0] for row in emissions_reinf]

EPD_reinf = cursor.execute("SELECT A1toA3_GWP FROM products "
                      "WHERE 'product_name [string]' LIKE '%Betonstahl%' "
                      "AND A1toA3_GWP IS NOT NULL "
                      "AND 'source [string]' NOT LIKE '%KBOB%' "
                      ).fetchall()
EPD_reinf_values = [row[0] for row in EPD_reinf]


print(EPD_reinf_values)

KBOB_reinf = cursor.execute("SELECT A1toA3_GWP FROM products "
                               "WHERE 'product_name [string]' LIKE '%Betonstahl%' "
                               "AND 'source [string]' LIKE '%KBOB%' "
                               ).fetchall()
KBOB_reinf_values = [row[0] for row in KBOB_reinf]

#------------------------------------------------------------------------------------------------------------------------
#plot concrete

# Combine all datasets to determine the range for bins
all_data = (emissions_concrete_values+EPD_concrete_values+ Ecoinvent_concrete_values + Betonsortenrechner_concrete_values)

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data, bins='auto')

# Create the filled histogram for emissions
plt.hist(emissions_concrete_values, bins=bins, label='Total', fill= False)

# Create the stacked histograms for Ecoinvent and Betonsortenrechner

plt.hist(Ecoinvent_concrete_values, bins=bins, stacked=True, label='Ecoinvent', color='gold', alpha=0.7)
plt.hist(Betonsortenrechner_concrete_values, bins=bins, stacked=True, label='Betonsortenrechner', color='lightseagreen', alpha=0.7)
plt.hist(EPD_concrete_values, bins=bins, stacked=True, alpha=0.7, label='Emissions', color='indianred', edgecolor='black')

#plt.hist(Ecoinvent, bins=3)
#plt.hist(Betonsortenrechner, bins=3)
plt.xlabel('A1-A3 GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Beton')
plt.legend()
#plt.axvline(KBOB_concrete_values, color = 'r')
plt.text(90, .1, 'KBOB Hochbaubeton', rotation=90, color = 'r')
plt.show()

#plot basic histogram
plt.boxplot([EPD_timber_values, KBOB_timber_values])
plt.show()

#------------------------------------------------------------------------------------------------------------------------
#plot wood

# Combine all datasets to determine the range for bins
all_data_timber = (emissions_timber_values+EPD_timber_values)

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_timber, bins='auto')

# Create the filled histogram for emissions
plt.hist(emissions_timber_values, bins=bins, label='Total', fill= False)

# Create the stacked histograms for Ecoinvent and Betonsortenrechner
plt.hist(EPD_timber_values, bins=bins, stacked=True, alpha=0.7, label='Emissions',  color='tan', edgecolor='black')

#plt.hist(Ecoinvent, bins=3)
#plt.hist(Betonsortenrechner, bins=3)
plt.xlabel('Total GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Holz')
plt.legend()
#plt.axvline(KBOB_timber_values[0], color = 'r')
#plt.text(255, 3, 'KBOB BSH CH', rotation=90, color = 'r')
#plt.axvline(KBOB_timber_values[1], color = 'r')
#plt.text(345, 3, 'KBOB BSH', rotation=90, color = 'r')
plt.show()


#------------------------------------------------------------------------------------------------------------------------
#plot Betonstahl

# Combine all datasets to determine the range for bins
all_data_reinf = (emissions_reinf_values+EPD_reinf_values)

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_reinf, bins='auto')

# Create the filled histogram for emissions
plt.hist(emissions_reinf_values, bins=bins, label='Total', fill= False)

# Create the stacked histograms for Ecoinvent and Betonsortenrechner
plt.hist(EPD_reinf_values, bins=bins, stacked=True, alpha=0.7, label='Emissions', color='lightsteelblue', edgecolor='black')


#plt.hist(Ecoinvent, bins=3)
#plt.hist(Betonsortenrechner, bins=3)
plt.xlabel('A1-A3 GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Betonstahl')
plt.legend()
#plt.axvline(KBOB_reinf_values, color = 'r')
#plt.text(775, 3, 'KBOB Betonstahl', rotation=90, color = 'r')
plt.show()



# #---------------
# emissions_reinf = cursor.execute(
#     "SELECT 'source [string]', A1toA3_GWP FROM products "
#     "WHERE 'product_name [string]' LIKE '%Betonstahl%' "
#     "AND A1toA3_GWP IS NOT NULL "
#     "AND A1toA3_GWP != 0"
# ).fetchall()
#
# # Creating a list with each entry as (source, A1toA3_GWP)
# emissions_list = [(entry[0], entry[1]) for entry in emissions_reinf]
#
# # Example print statement to visualize the output
# print(emissions_list)
#
#
# # Unpacking the emissions_list into two separate lists
# sources, gwps = zip(*emissions_list)
#
# # Create a horizontal bar plot
# plt.barh(sources, gwps, color='lightsteelblue')
#
# # Add labels and title
# plt.xlabel('A1-A3 GWP [kg CO2-eq/t]')
# plt.ylabel("source")
# plt.title("Betonstahl")
#
# # Show the plot
# plt.tight_layout()
# plt.show()