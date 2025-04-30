import matplotlib.pyplot as plt
import sqlite3
import numpy as np
import statistics

# create or open database sustainability
connection = sqlite3.connect('database_250326.db')
# create cursor object
cursor = connection.cursor()
#------------------------------------------------------------------------------------------------------------------------
#extract values for concrete
#
emissions_concrete = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    "WHERE type LIKE '%ready mixed concrete%' "
                    "AND A1toA3_GWP IS NOT NULL "
                    "AND A1toA3_GWP != 0 "
                    "AND type LIKE '%ready mixed concrete%' "
                    ).fetchall()
emissions_concrete_values = [row[0] for row in emissions_concrete]

EPD_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%Betonsortenrechner%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND "source [string]" NOT LIKE '%KBOB%'
                    """).fetchall()
EPD_concrete_values = [row[0] for row in EPD_concrete]

KBOB_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND "source [string]" LIKE '%KBOB%'
                    """).fetchall()
KBOB_concrete_values = [row[0] for row in KBOB_concrete]


Ecoinvent_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND "source [string]" LIKE '%Ecoinvent%'
                    AND A1toA3_GWP != 0
                    """).fetchall()
Ecoinvent_concrete_values = [row[0] for row in Ecoinvent_concrete]

Betonsortenrechner_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND "source [string]" LIKE '%Betonsortenrechner%'
                    """).fetchall()
Betonsortenrechner_concrete_values = [row[0] for row in Betonsortenrechner_concrete]

print(EPD_concrete_values)
print(np.quantile(EPD_concrete_values, 0.1))
#------------------------------------------------------------------------------------------------------------------------
#extract values for wood
emissions_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS  NOT NULL
                    AND "source [string]" NOT LIKE '%Studiengemeinschaft%'
                    """).fetchall()
emissions_timber_values = [row[0] for row in emissions_timber]

EPD_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%KBOB%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND EPD_ID NOT LIKE '%verifizierung%'
                    """).fetchall()
EPD_timber_values = [row[0] for row in EPD_timber]


Ecoinvent_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS NOT NULL
                    AND "source [string]" LIKE '%Ecoinvent%'
                    """).fetchall()
Ecoinvent_timber_values = [row[0] for row in Ecoinvent_timber]


KBOB_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND "source [string]" LIKE '%KBOB%'
                    """).fetchall()
KBOB_timber_values = [row[0] for row in KBOB_timber]

#------------------------------------------------------------------------------------------------------------------------
#extract values for reinforcement

emissions_reinf = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE type LIKE '%steel reinforcing bar%'
                        AND Total_GWP IS NOT NULL
                        AND A1toA3_GWP != 0
                        """).fetchall()
emissions_reinf_values = [row[0] for row in emissions_reinf]

EPD_reinf = cursor.execute("""
                      SELECT Total_GWP FROM products
                      WHERE type LIKE '%steel reinforcing bar%'
                      AND Total_GWP IS NOT NULL
                      AND "source [string]" NOT LIKE '%KBOB%'
                      """).fetchall()
EPD_reinf_values = [row[0] for row in EPD_reinf]

KBOB_reinf = cursor.execute(
                            """
                            SELECT Total_GWP FROM products
                            WHERE type LIKE '%steel reinforcing bar%'
                            AND "source [string]" LIKE '%KBOB%'
                            """).fetchall()
KBOB_reinf_values = [row[0] for row in KBOB_reinf]

reinf_min = min(EPD_reinf_values)
reinf_max =max(EPD_reinf_values)

#------------------------------------------------------------------------------------------------------------------------
#extract values for steel
emissions_steel = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE type LIKE '%structural steel profile%'
                        AND Total_GWP IS NOT NULL
                        AND Total_GWP != 0
                        """).fetchall()
emissions_steel_values = [row[0] for row in emissions_steel]

EPD_steel = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE type LIKE '%structural steel profile%'
                        AND Total_GWP IS NOT NULL
                        AND "source [string]" NOT LIKE '%KBOB%'
                        """).fetchall()
EPD_steel_values = [row[0] for row in EPD_reinf]

KBOB_steel = cursor.execute("""
                        SELECT A1toA3_GWP FROM products
                        WHERE type LIKE '%structural steel profile%'
                        AND "source [string]" LIKE '%KBOB%'
                        """).fetchall()
KBOB_steel_values = [row[0] for row in KBOB_steel]

#------------------------------------------------------------------------------------------------------------------------
#plot concrete
import numpy as np
import matplotlib.pyplot as plt

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))#, sharey=True)

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
ax1.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='forestgreen', edgecolor='black')
ax1.scatter(bin_centers_Ecoinvent, hist_Ecoinvent, alpha=0.3, label='Ecoinvent', color='gray', edgecolor='black')
ax1.scatter(bin_centers_Betonsortenrechner, hist_Betonsortenrechner, alpha=0.3, label='Betonsortenrechner', color='plum', edgecolor='black')

# Set y-axis to display only integer values
ymin = 0
ymax = 11
ax1.set_yticks(range(ymin, ymax))
ax1.set_ylim(ymin,ymax)


ax1.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax1.set_ylabel('#')
ax1.set_title('Beton')
ax1.legend(loc='upper right')

# Add vertical lines and text for KBOB values
ax1.axvline(KBOB_concrete_values[0], linestyle='--', alpha=0.5, color='tomato')
ax1.text(KBOB_concrete_values[0]+0.5, max(hist_EPD)/2, 'KBOB Hochbaubeton',alpha=0.5, rotation=90, color='tomato')

ax1.axvline(statistics.mean(EPD_concrete_values), linestyle='--', alpha=0.6, color='forestgreen')
ax1.text(statistics.mean(EPD_concrete_values)+0.5, ymin+0.5, r'$\mu$', color = 'forestgreen', alpha = 0.6)
ax1.fill_betweenx([ymin, ymax], np.quantile(EPD_concrete_values,0.1), np.quantile(EPD_concrete_values,0.9), color='forestgreen', alpha=0.05)
ax1.text(np.quantile(EPD_concrete_values, 0.1)+2, ymax-0.5, '90 %', color = 'forestgreen', alpha=0.5)

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
print(hist_EPD)

# Create the scatter plots
colors = ['white' if epd < 0 else 'peru' for epd in bin_centers_EPD]
ax2.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions',
            c=colors, edgecolor='black')

#ax2.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='peru', edgecolor='black')
ax2.scatter(bin_centers_Ecoinvent, hist_Ecoinvent, alpha=0.3, label='Ecoinvent Emissions', color='gray', edgecolor='black')

ax2.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax2.set_ylabel('#')
ax2.set_title('Holz')
ax2.legend(loc='upper right')

ax2.set_yticks(range(ymin, ymax))
ax2.set_ylim(ymin,ymax)

ax2.set_xlim(-300)

# Add vertical lines and text for KBOB values
ax2.axvline(KBOB_timber_values[0], linestyle='--', alpha= 0.5, color='tomato')
ax2.text(KBOB_timber_values[0]+15, max(max(hist_Ecoinvent), max(hist_EPD))/2, 'KBOB BSH CH', rotation=90, alpha=0.5, color='tomato')
ax2.axvline(KBOB_timber_values[1], linestyle='--', alpha=0.5, color='coral')
ax2.text(KBOB_timber_values[1]+15, max(max(hist_Ecoinvent), max(hist_EPD))/2, 'KBOB BSH', rotation=90, alpha = 0.5, color='darkorange')

EPD_timber_values_pos = [x for x in EPD_timber_values if x >= 0]

ax2.axvline(statistics.mean(EPD_timber_values_pos), linestyle='--', alpha=0.6, color='peru')
ax2.text(statistics.mean(EPD_timber_values_pos)+15, ymin+0.5, r'$\mu$', color = 'peru', alpha = 0.6)
ax2.fill_betweenx([ymin, ymax], np.quantile(EPD_timber_values_pos,0.1), np.quantile(EPD_timber_values_pos,0.9), color='peru', alpha=0.1)
ax2.text(np.quantile(EPD_timber_values_pos, 0.1)+20, ymax-0.5, '90 %', color = 'peru', alpha=0.5)

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

colors = ['red' if epd < 0 else 'steelblue' for epd in hist_EPD]

ax3.scatter(bin_centers, hist_EPD, alpha=0.7, label='EPD Emissions',
            color=colors, edgecolor='black')

# ax3.scatter(bin_centers, hist_EPD_pos, alpha=0.7, label='EPD Emissions', color='steelblue', edgecolor='black')
# ax3.scatter(bin_centers, hist_EPD_neg, alpha=0.2, label='EPD Emissions', color='steelblue', edgecolor='black')

ax3.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax3.set_ylabel('#')
ax3.set_title('Betonstahl')
ax3.legend(loc='upper right')

ax3.set_yticks(range(ymin, ymax))
ax3.set_ylim(ymin,ymax)

# Add vertical lines and text for KBOB values
ax3.axvline(KBOB_reinf_values[0], linestyle='--', alpha = 0.5, color='tomato')
ax3.text(KBOB_reinf_values[0]+0.5, max(hist_EPD)/2, 'KBOB Bewehrung', rotation=90, alpha= 0.5, color='tomato')

ax3.axvline(368, color='blue')
ax3.text(368+0.5, max(hist_EPD)/2, 'Stahl Gerlafingen', rotation=90, color='blue')

# Add vertical lines and text for KBOB values
ax3.axvline(statistics.mean(EPD_reinf_values), linestyle='--', alpha=0.6, color='blue')
ax3.text(statistics.mean(EPD_reinf_values)+15, ymin+0.5, r'$\mu$', color = 'blue', alpha = 0.6)
ax3.fill_betweenx([ymin, ymax], np.quantile(EPD_reinf_values,0.1), np.quantile(EPD_reinf_values,0.9), color='blue', alpha=0.1)
ax3.text(np.quantile(EPD_reinf_values, 0.1)+10, ymax-0.5, '90 %', color = 'blue', alpha=0.5)

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
ax4.scatter(bin_centers, hist_EPD, alpha=0.7, label='EPD Emissions', color='deepskyblue', edgecolor='black')

ax4.set_yticks(range(ymin, ymax))
ax4.set_ylim(ymin,ymax)


ax4.set_xlabel('Total GWP [kg CO2-eq/t]')
ax4.set_ylabel('#')
ax4.set_title('Baustahl')
ax4.legend(loc='upper right')



# Add vertical lines and text for KBOB values
ax4.axvline(KBOB_steel_values[0], linestyle='--', alpha = 0.5, color='tomato')
ax4.text(KBOB_steel_values[0]+10, max(hist_EPD)/2, 'KBOB Baustahl', rotation=90, alpha=0.5, color='tomato')

# Add vertical lines and text for KBOB values
ax4.axvline(statistics.mean(EPD_steel_values), linestyle='--', alpha=0.6, color='deepskyblue')
ax4.text(statistics.mean(EPD_steel_values)+15, ymin+0.5, r'$\mu$', color = 'deepskyblue', alpha = 0.6)
ax4.fill_betweenx([ymin, ymax], np.quantile(EPD_steel_values,0.1), np.quantile(EPD_steel_values,0.9), color='deepskyblue', alpha=0.1)
ax4.text(np.quantile(EPD_steel_values, 0.1)+10, ymax-0.5, '90 %', color = 'deepskyblue', alpha=0.5)

plt.show()

#--------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import statistics


# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 5), sharey=True)

ax1.boxplot([EPD_reinf_values])

# Adding scatter points for the data
y = EPD_concrete_values
x = np.random.normal(1, 0.04, len(y))
ax1.plot(x, y, 'r.', alpha=0.5)

# Adding titles and labels
ax1.set_xlabel('Gruppe')
ax1.set_ylabel('Total GWP [kg CO$_2$-eq/t]')
ax1.set_title('Beton')

# Customizing x-axis labels
ax1.set_xticks([1], ['EPD Emissions'])

ax2.boxplot([EPD_reinf_values])
# Adding scatter points for the data
y = EPD_timber_values
x = np.random.normal(1, 0.04, len(y))
ax2.plot(x, y, 'r.', alpha=0.5)

# Adding titles and labels
ax2.set_xlabel('Gruppe')
ax2.set_ylabel('Total GWP [kg CO2-eq/t]')
ax2.set_title('Holz')

# Customizing x-axis labels
ax2.set_xticks([1], ['EPD Emissions'])


ax3.boxplot([EPD_reinf_values])
ax4.boxplot([EPD_reinf_values])
# Display the plot
plt.show()
