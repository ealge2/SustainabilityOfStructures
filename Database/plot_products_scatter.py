import matplotlib.pyplot as plt
import sqlite3
import numpy as np
import statistics

# create or open database sustainability
connection = sqlite3.connect('data.db')
# create cursor object
cursor = connection.cursor()
#------------------------------------------------------------------------------------------------------------------------
#extract values for concrete
emissions_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material LIKE '%concrete%' "
                    "AND A1toA3_GWP IS NOT NULL "
                    "AND A1toA3_GWP != 0 "
                    ).fetchall()
emissions_concrete_values = [row[0] for row in emissions_concrete]

EPD_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material LIKE '%concrete%' "
                    "AND A1toA3_GWP IS NOT NULL "
                    "AND source NOT LIKE '%Betonsortenrechner%' "
                    "AND source NOT LIKE '%Ecoinvent%' "
                    "AND source NOT LIKE '%KBOB%' "
                    ).fetchall()
EPD_concrete_values = [row[0] for row in EPD_concrete]

KBOB_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material LIKE '%concrete%' "
                    "AND source LIKE '%KBOB%' "
                    ).fetchall()
KBOB_concrete_values = [row[0] for row in KBOB_concrete]

Ecoinvent_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material LIKE '%concrete%' "
                    "AND source LIKE '%Ecoinvent%' "
                    "AND A1toA3_GWP != 0 "
                    ).fetchall()
Ecoinvent_concrete_values = [row[0] for row in Ecoinvent_concrete]

Betonsortenrechner_concrete = cursor.execute(
                    "SELECT A1toA3_GWP FROM products "
                    "WHERE material LIKE '%concrete%' "
                    "AND source LIKE '%Betonsortenrechner%' "
                    ).fetchall()
Betonsortenrechner_concrete_values = [row[0] for row in Betonsortenrechner_concrete]

#------------------------------------------------------------------------------------------------------------------------
#extract values for wood
emissions_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    "WHERE material LIKE '%timber%' "
                    "AND Total_GWP IS  NOT NULL "
                    "AND source NOT LIKE '%Studiengemeinschaft%' "
                    ).fetchall()
emissions_timber_values = [row[0] for row in emissions_timber]

EPD_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    " WHERE material LIKE '%timber%' "
                    "AND Total_GWP IS NOT NULL "
                    "AND source NOT LIKE '%KBOB%' "
                    "AND source NOT LIKE '%Ecoinvent%' "
                    "AND ï»¿EPD_ID NOT LIKE '%verifizierung%' "
                    ).fetchall()
EPD_timber_values = [row[0] for row in EPD_timber]


Ecoinvent_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    " WHERE material LIKE '%timber%' "
                    "AND Total_GWP IS NOT NULL "
                    "AND source LIKE '%Ecoinvent%' "
                    ).fetchall()
Ecoinvent_timber_values = [row[0] for row in Ecoinvent_timber]


KBOB_timber = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    "WHERE material LIKE '%timber%' "
                    "AND source LIKE '%KBOB%' "
                    ).fetchall()
KBOB_timber_values = [row[0] for row in KBOB_timber]

#------------------------------------------------------------------------------------------------------------------------
#extract values for reinforcement

emissions_reinf = cursor.execute("SELECT A1toA3_GWP FROM products "
                           " WHERE product_name LIKE '%Betonstahl%' "
                           "AND A1toA3_GWP IS NOT NULL "
                           "AND A1toA3_GWP != 0 "
                           ).fetchall()
emissions_reinf_values = [row[0] for row in emissions_reinf]

EPD_reinf = cursor.execute("SELECT A1toA3_GWP FROM products "
                      "WHERE product_name LIKE '%Betonstahl%' "
                      "AND A1toA3_GWP IS NOT NULL "
                      "AND source NOT LIKE '%KBOB%' "
                      ).fetchall()
EPD_reinf_values = [row[0] for row in EPD_reinf]

KBOB_reinf = cursor.execute("SELECT A1toA3_GWP FROM products "
                               "WHERE product_name LIKE '%Betonstahl%' "
                               "AND source LIKE '%KBOB%' "
                               ).fetchall()
KBOB_reinf_values = [row[0] for row in KBOB_reinf]

#------------------------------------------------------------------------------------------------------------------------
#extract values for steel
emissions_steel = cursor.execute("SELECT A1toA3_GWP FROM products "
                           " WHERE product_name LIKE '%Baustahl%' "
                           "AND A1toA3_GWP IS NOT NULL "
                           "AND A1toA3_GWP != 0 "
                           ).fetchall()
emissions_steel_values = [row[0] for row in emissions_steel]

EPD_steel = cursor.execute("SELECT A1toA3_GWP FROM products "
                      "WHERE product_name LIKE '%Baustahl%' "
                      "AND A1toA3_GWP IS NOT NULL "
                      "AND source NOT LIKE '%KBOB%' "
                      ).fetchall()
EPD_steel_values = [row[0] for row in EPD_reinf]

KBOB_steel = cursor.execute("SELECT A1toA3_GWP FROM products "
                               "WHERE product_name LIKE '%Baustahl%' "
                               "AND source LIKE '%KBOB%' "
                               ).fetchall()
KBOB_steel_values = [row[0] for row in KBOB_steel]

#------------------------------------------------------------------------------------------------------------------------
#plot concrete

import numpy as np
import matplotlib.pyplot as plt

# Combine all datasets to determine the range for bins
all_data_concrete = Ecoinvent_concrete_values + Betonsortenrechner_concrete_values + EPD_concrete_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_concrete, bins='auto')

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histograms
hist_Ecoinvent, _ = np.histogram(Ecoinvent_concrete_values, bins=bins)
hist_Betonsortenrechner, _ = np.histogram(Betonsortenrechner_concrete_values, bins=bins)
hist_EPD, _ = np.histogram(EPD_concrete_values, bins=bins)

# Filter out zero frequency bins
non_zero_Ecoinvent = hist_Ecoinvent > 0
non_zero_Betonsortenrechner = hist_Betonsortenrechner > 0
non_zero_EPD = hist_EPD > 0

bin_centers_Ecoinvent = bin_centers[non_zero_Ecoinvent]
bin_centers_Betonsortenrechner = bin_centers[non_zero_Betonsortenrechner]
bin_centers_EPD = bin_centers[non_zero_EPD]

hist_Ecoinvent = hist_Ecoinvent[non_zero_Ecoinvent]
hist_Betonsortenrechner = hist_Betonsortenrechner[non_zero_Betonsortenrechner]
hist_EPD = hist_EPD[non_zero_EPD]

# Create the scatter plots
plt.scatter(bin_centers_Ecoinvent, hist_Ecoinvent, alpha=0.7, label='Ecoinvent', color='gray', edgecolor='black')
plt.scatter(bin_centers_Betonsortenrechner, hist_Betonsortenrechner, alpha=0.7, label='Betonsortenrechner', color='plum', edgecolor='black')
plt.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='lightseagreen', edgecolor='black')

plt.xlabel('Total GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Beton')
plt.legend(loc='upper right')


# Add vertical lines and text for KBOB values
plt.axvline(KBOB_concrete_values[0], color='tomato')
plt.text(KBOB_concrete_values[0]+0.5, max(max(hist_Ecoinvent), max(hist_Betonsortenrechner), max(hist_EPD))/2, 'KBOB Hochbaubeton', rotation=90, color='tomato')



plt.show()

'''# Plot a boxplot for EPD_timber_values and KBOB_timber_values
plt.boxplot([EPD_timber_values, KBOB_timber_values])
plt.show()'''

#------------------------------------------------------------------------------------------------------------------------
#plot wood

# Combine all datasets to determine the range for bins
all_data_timber = EPD_timber_values + Ecoinvent_timber_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_timber, bins='auto')

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histograms
hist_EPD, _ = np.histogram(EPD_timber_values, bins=bins)
hist_Ecoinvent, _ = np.histogram(Ecoinvent_timber_values, bins=bins)

# Filter out zero frequency bins
non_zero_Ecoinvent = hist_Ecoinvent > 0
non_zero_EPD = hist_EPD > 0

bin_centers_Ecoinvent = bin_centers[non_zero_Ecoinvent]
bin_centers_EPD = bin_centers[non_zero_EPD]

hist_Ecoinvent = hist_Ecoinvent[non_zero_Ecoinvent]
hist_EPD = hist_EPD[non_zero_EPD]


# Create the scatter plots
plt.scatter(bin_centers_Ecoinvent, hist_Ecoinvent, label='Ecoinvent Emissions', color='gray', edgecolor='black')
plt.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='gold', edgecolor='black')

plt.xlabel('Total GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Holz')
plt.legend(loc='upper right')

# Set y-axis to display only integer values
plt.yticks(range(int(np.floor(min(hist_EPD))), int(np.ceil(max(hist_EPD)))+1))

# Add vertical lines and text for KBOB values
plt.axvline(KBOB_timber_values[0], color='tomato')
plt.text(KBOB_timber_values[0]+0.5, max(max(hist_Ecoinvent), max(hist_EPD))/2, 'KBOB BSH CH', rotation=90, color='tomato')
plt.axvline(KBOB_timber_values[1], color='coral')
plt.text(KBOB_timber_values[1]+0.5, max(max(hist_Ecoinvent), max(hist_EPD))/2, 'KBOB BSH', rotation=90, color='darkorange')

plt.show()

#------------------------------------------------------------------------------------------------------------------------
#Betonstahl

# Combine all datasets to determine the range for bins
all_data_reinf = EPD_reinf_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_reinf, bins='auto')

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histogram for EPD_timber_values
hist_EPD, _ = np.histogram(EPD_reinf_values, bins=bins)


# Create the scatter plot for EPD values
plt.scatter(bin_centers, hist_EPD, alpha=0.7, label='EPD Emissions', color='steelblue', edgecolor='black')

plt.xlabel('Total GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Betonstahl')
plt.legend(loc='upper right')

# Set y-axis to display only integer values
plt.yticks(range(int(np.floor(min(hist_EPD))), int(np.ceil(max(hist_EPD)))+1))

# Add vertical lines and text for KBOB values
plt.axvline(KBOB_reinf_values[0], color='tomato')
plt.text(KBOB_reinf_values[0]+0.5, max(hist_EPD)/2, 'KBOB Bewehrung', rotation=90, color='tomato')

plt.axvline(368, color='blue')
plt.text(368+0.5, max(hist_EPD)/2, 'Stahl Gerlafingen', rotation=90, color='blue')

plt.show()
#------------------------------------------------------------------------------------------------------------------------
#plot steel
# Combine all datasets to determine the range for bins
all_data_steel = EPD_steel_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_steel, bins='auto')

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histogram for EPD_timber_values
hist_EPD, _ = np.histogram(EPD_steel_values, bins=bins)

# Create the scatter plot for EPD values
plt.scatter(bin_centers, hist_EPD, alpha=0.7, label='EPD Emissions', color='cornflowerblue', edgecolor='black')

plt.xlabel('Total GWP [kg CO2-eq/t]')
plt.ylabel('#')
plt.title('Baustahl')
plt.legend(loc='upper right')

# Set y-axis to display only integer values
plt.yticks(range(int(np.floor(min(hist_EPD))), int(np.ceil(max(hist_EPD)))+1))

# Add vertical lines and text for KBOB values
plt.axvline(KBOB_steel_values[0], color='tomato')
plt.text(KBOB_steel_values[0]+0.5, max(hist_EPD)/2, 'KBOB Baustahl', rotation=90, color='tomato')

plt.show()