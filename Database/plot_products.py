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
EPD_steel_values = [row[0] for row in EPD_steel]

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


plt.rcParams.update({'font.size': 16})


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
#ax1.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='forestgreen', edgecolor='black')
#ax1.scatter(bin_centers_Ecoinvent, hist_Ecoinvent, alpha=0.3, label='Ecoinvent', color='gray', edgecolor='black')
#ax1.scatter(bin_centers_Betonsortenrechner, hist_Betonsortenrechner, alpha=0.3, label='Betonsortenrechner', color='plum', edgecolor='black')
bin_widths = np.diff(bins)
ax1.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color='forestgreen', edgecolor='black', alpha=0.4)
ax1.scatter(EPD_concrete_values, np.random.normal(2, 0.5, len(EPD_concrete_values)), alpha=0.7, label='EPD Emissions', color = 'white', edgecolor='black')

ymin = 0
ymax = 14
ax1.set_yticks(range(ymin, ymax))
ax1.set_ylim(ymin,ymax)

ax1.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax1.set_ylabel('#')
ax1.set_title('Beton')
#ax1.legend(loc='upper right')

# Add vertical lines and text for KBOB values
ax1.axvline(KBOB_concrete_values[0], linestyle='--', alpha=0.9, color='tomato')
ax1.text(KBOB_concrete_values[0]+0.2, 6, 'KBOB Hochbaubeton', fontsize = 12, alpha=0.9, rotation=90, color='tomato')

ax1.axvline(statistics.mean(EPD_concrete_values), linestyle='--', alpha=0.9, color='forestgreen')
ax1.text(statistics.mean(EPD_concrete_values), -0.7, r'$\mu$', color = 'forestgreen', alpha = 0.9)
ax1.fill_betweenx([ymin, ymax], np.quantile(EPD_concrete_values,0.1), np.quantile(EPD_concrete_values,0.9), color='grey', alpha=0.1, zorder=1)
ax1.text(np.quantile(EPD_concrete_values, 0.1)+1, ymax-1, '80%', color = 'forestgreen', alpha=0.9)

#------------------------------------------------------------------------------------------------------------------------
#plot wood

# Combine all datasets to determine the range for bins
all_data_timber = EPD_timber_values #+ Ecoinvent_timber_values

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
colors = ['white' if epd < 0 else 'peru' for epd in bin_centers_EPD]
bin_widths = np.diff(bins)
ax2.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color=colors, edgecolor='black', alpha=0.4)
ax2.scatter(EPD_timber_values, np.random.normal(2, 0.5, len(EPD_timber_values)), alpha=0.7, label='EPD Emissions', color = 'white', edgecolor='black')

#ax2.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='peru', edgecolor='black')
#ax2.scatter(bin_centers_Ecoinvent, hist_Ecoinvent, alpha=0.3, label='Ecoinvent Emissions', color='gray', edgecolor='black')

ax2.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax2.set_ylabel('#')
ax2.set_title('Holz')
#ax2.legend(loc='upper right')

ax2.set_yticks(range(ymin, ymax))
ax2.set_ylim(ymin,ymax)

ax2.set_xlim(-300)

# Add vertical lines and text for KBOB values
ax2.axvline(KBOB_timber_values[0], linestyle='--', alpha= 0.9, color='tomato')
ax2.text(KBOB_timber_values[0]+10, 8, 'KBOB BSH CH', fontsize = 12, rotation=90, alpha=0.9, color='tomato')
ax2.axvline(KBOB_timber_values[1], linestyle='--', alpha=0.9, color='darkorange')
ax2.text(KBOB_timber_values[1]+10, 8, 'KBOB BSH', fontsize = 12, rotation=90, alpha = 0.9, color='darkorange')
ax2.axvline(288, linestyle='--', alpha=0.9, color='coral')
ax2.text(288+10, 8, 'KBOB KVH', fontsize = 12, rotation=90, alpha = 0.9, color='coral')

EPD_timber_values_pos = [x for x in EPD_timber_values if x >= 0]

ax2.axvline(statistics.mean(EPD_timber_values_pos), linestyle='--', alpha=0.9, color='peru')
ax2.text(statistics.mean(EPD_timber_values_pos), -0.7, r'$\mu$', color = 'peru', alpha = 0.9)
ax2.fill_betweenx([ymin, ymax], np.quantile(EPD_timber_values_pos,0.1), np.quantile(EPD_timber_values_pos,0.9), color='grey', alpha=0.1, zorder=1)
ax2.text(np.quantile(EPD_timber_values_pos, 0.1)-3, ymax-1, '80%', color = 'peru', alpha=0.9)

#------------------------------------------------------------------------------------------------------------------------
#Betonstahl
# Combine all datasets to determine the range for bins
all_data_reinf = EPD_reinf_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_reinf, bins='auto')

# Calculate the bin centers
bin_centers_EPD = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histogram for EPD_timber_values
hist_EPD, _ = np.histogram(EPD_reinf_values, bins=bins)

# Create the scatter plot for EPD values

colors = ['red' if epd < 0 else 'steelblue' for epd in hist_EPD]

bin_widths = np.diff(bins)
ax3.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color='blue', edgecolor='black', alpha=0.3)
#ax3.scatter(bin_centers_EPD, hist_EPD, alpha=0.7, label='EPD Emissions', color='blue', edgecolor='black')
ax3.scatter(EPD_reinf_values, np.random.normal(2, 0.5, len(EPD_reinf_values)), alpha=0.7, label='EPD Emissions', color = 'white', edgecolor='black')
# ax3.scatter(bin_centers, hist_EPD_pos, alpha=0.7, label='EPD Emissions', color='steelblue', edgecolor='black')
# ax3.scatter(bin_centers, hist_EPD_neg, alpha=0.2, label='EPD Emissions', color='steelblue', edgecolor='black')

ax3.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax3.set_ylabel('#')
ax3.set_title('Betonstahl')
#ax3.legend(loc='upper right')

ax3.set_yticks(range(ymin, ymax))
ax3.set_ylim(ymin,ymax)

# Add vertical lines and text for KBOB values
ax3.axvline(KBOB_reinf_values[0], linestyle='--', alpha = 0.7, color='tomato')
ax3.text(KBOB_reinf_values[0]+5, 6, 'KBOB Bewehrung', fontsize = 12, rotation=90, alpha= 0.9, color='tomato')

ax3.axvline(368, color='grey',linestyle='dotted')
ax3.text(368+10, 4, 'Stahl Gerlafingen', fontsize = 12, rotation=90, color='grey')

# Add vertical lines and text for KBOB values
ax3.axvline(statistics.mean(EPD_reinf_values), linestyle='--', alpha=0.9, color='blue')
ax3.text(statistics.mean(EPD_reinf_values)-10, -0.7, r'$\mu$', color = 'blue', alpha = 0.9)
ax3.fill_betweenx([ymin, ymax], np.quantile(EPD_reinf_values,0.1), np.quantile(EPD_reinf_values,0.9), color='grey', alpha=0.1, zorder=1)
ax3.text(np.quantile(EPD_reinf_values, 0.1)+10, ymax-1, '80%', color = 'blue', alpha=0.9)

#------------------------------------------------------------------------------------------------------------------------
#plot steel
# Combine all datasets to determine the range for bins
all_data_steel = EPD_steel_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_steel, bins='auto')

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])
# Calculate the bin centers
bin_centers_EPD = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histogram for EPD_timber_values
hist_EPD, _ = np.histogram(EPD_steel_values, bins=bins)

# Create the scatter plot for EPD values
#ax4.scatter(bin_centers, hist_EPD, alpha=0.7, label='EPD Emissions', color='lightskyblue', edgecolor='black')
bin_widths = np.diff(bins)
ax4.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color='lightskyblue', edgecolor='black', alpha=0.4)
ax4.scatter(EPD_steel_values, np.random.normal(2, 0.5, len(EPD_steel_values)), alpha=0.7, label='EPD Emissions', color = 'white', edgecolor='black')

ax4.set_yticks(range(ymin, ymax))
ax4.set_ylim(ymin,ymax)


ax4.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax4.set_ylabel('#')
ax4.set_title('Baustahl')
#ax4.legend(loc='upper right')



# Add vertical lines and text for KBOB values
ax4.axvline(KBOB_steel_values[0], linestyle='--', alpha = 0.7, color='tomato')
ax4.text(KBOB_steel_values[0]+20, 6, 'KBOB Baustahl', fontsize = 12, rotation=90, alpha=0.9, color='tomato')

# Add vertical lines and text for KBOB values
ax4.axvline(statistics.mean(EPD_steel_values), linestyle='--', alpha=0.9, color='lightskyblue')
ax4.text(statistics.mean(EPD_steel_values), -0.7, r'$\mu$', color = 'lightskyblue', alpha = 0.9)
ax4.fill_betweenx([ymin, ymax], np.quantile(EPD_steel_values,0.1), np.quantile(EPD_steel_values,0.9), color='grey', alpha=0.1, zorder=1)
ax4.text(np.quantile(EPD_steel_values, 0.1)+10, ymax-1, '80%', color = 'lightskyblue', alpha=0.9)

plt.show()

#--------------------------------------------
# import matplotlib.pyplot as plt
# import numpy as np
# import statistics
#
# # Create subplots
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 5))
#
# ax1.boxplot([EPD_concrete_values])# vert=False)
#
# # Adding scatter points for the data
# y = EPD_concrete_values
# x = np.random.normal(1, 0.04, len(y))
# ax1.plot(x, y, 'r.', color = 'forestgreen', alpha=0.5)
#
# # Adding titles and labels
# ax1.set_xlabel('Beton')
# ax1.set_ylabel('Total GWP [kg CO$_2$-eq/t]')
# ax1.set_title('Beton')
#
# # Customizing x-axis labels
# ax1.set_xticks([1], ['EPD Emissions'])
#
# ax2.boxplot([EPD_timber_values])# vert=False)
# # Adding scatter points for the data
# y = EPD_timber_values
# x = np.random.normal(1, 0.04, len(y))
# ax2.plot(x, y, 'r.', color = 'peru', alpha=0.5)
#
# # Adding titles and labels
# ax2.set_xlabel('Gruppe')
# ax2.set_ylabel('Total GWP [kg CO$_2$-eq/t]')
# ax2.set_title('Holz')
#
# # Customizing x-axis labels
# ax2.set_xticks([1], ['EPD Emissions'])
#
#
# ax3.boxplot([EPD_reinf_values])
#
# # Adding scatter points for the data
# y = EPD_reinf_values
# x = np.random.normal(1, 0.04, len(y))
# ax3.plot(x, y, 'r.', color = 'blue', alpha=0.5)
#
# # Adding titles and labels
# ax3.set_xlabel('Gruppe')
# ax3.set_ylabel('Total GWP [kg CO$_2$-eq/t]')
# ax3.set_title('Bewehrungsstahl')
#
# # Customizing x-axis labels
# ax3.set_xticks([1], ['EPD Emissions'])
#
# ax4.boxplot([EPD_steel_values])
#
# # Adding scatter points for the data
# y = EPD_steel_values
# x = np.random.normal(1, 0.04, len(y))
# ax4.plot(x, y, 'r.', color = 'deepskyblue', alpha=0.5)
#
# # Adding titles and labels
# ax4.set_xlabel('Gruppe')
# ax4.set_ylabel('Total GWP [kg CO$_2$-eq/t]')
# ax4.set_title('Baustahl')
# # Display the plot
# plt.show()
