import matplotlib.pyplot as plt
import sqlite3
import numpy as np
import statistics

# create or open database sustainability
connection = sqlite3.connect('database_251030.db')
# create cursor object
cursor = connection.cursor()
#------------------------------------------------------------------------------------------------------------------------
#extract values for concrete
emissions_concrete = cursor.execute(
                    "SELECT Total_GWP FROM products "
                    "WHERE MATERIAL LIKE '%ready mixed concrete%' "
                    "AND A1toA3_GWP IS NOT NULL "
                    "AND A1toA3_GWP != 0 "
                    "AND MATERIAL LIKE '%ready mixed concrete%' "
                    ).fetchall()
emissions_concrete_values = [row[0] for row in emissions_concrete]

EPD_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL LIKE '%ready_mixed_concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                    AND "SOURCE" NOT LIKE '%Ecoinvent%'
                    AND "SOURCE" NOT LIKE '%KBOB%'
                    """).fetchall()
EPD_concrete_values = [row[0] for row in EPD_concrete]

KBOB_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL LIKE '%ready_mixed_concrete%'
                    AND "SOURCE" LIKE '%KBOB%'
                    """).fetchall()
KBOB_concrete_values = [row[0] for row in KBOB_concrete]

Ecoinvent_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL LIKE '%ready_mixed_concrete%'
                    AND "SOURCE" LIKE '%Ecoinvent%'
                    AND A1toA3_GWP != 0
                    """).fetchall()
Ecoinvent_concrete_values = [row[0] for row in Ecoinvent_concrete]

Betonsortenrechner_concrete = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL LIKE '%ready_mixed_concrete%'
                    AND "SOURCE" LIKE '%Betonsortenrechner%'
                    """).fetchall()
Betonsortenrechner_concrete_values = [row[0] for row in Betonsortenrechner_concrete]

print(EPD_concrete_values)
print(np.quantile(EPD_concrete_values, 0.1))
#------------------------------------------------------------------------------------------------------------------------
#extract values for wood
emissions_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL_GROUP LIKE '%wood%'
                    AND Total_GWP IS  NOT NULL
                    AND "SOURCE" NOT LIKE '%Studiengemeinschaft%'
                    """).fetchall()
emissions_timber_values = [row[0] for row in emissions_timber]

EPD_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL_GROUP LIKE '%wood%'
                    AND Total_GWP IS NOT NULL
                    AND "SOURCE" NOT LIKE '%KBOB%'
                    AND "SOURCE" NOT LIKE '%Ecoinvent%'
                    AND EPD_ID NOT LIKE '%verifizierung%'
                    """).fetchall()
EPD_timber_values = [row[0] for row in EPD_timber]


Ecoinvent_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL_GROUP LIKE '%wood%'
                    AND Total_GWP IS NOT NULL
                    AND "SOURCE" LIKE '%Ecoinvent%'
                    """).fetchall()
Ecoinvent_timber_values = [row[0] for row in Ecoinvent_timber]


KBOB_timber = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE MATERIAL_GROUP LIKE '%wood%'
                    AND "SOURCE" LIKE '%KBOB%'
                    """).fetchall()
KBOB_timber_values = [row[0] for row in KBOB_timber]

#------------------------------------------------------------------------------------------------------------------------
#extract values for reinforcement

emissions_reinf = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE MATERIAL LIKE '%steel_reinforcing_bar%'
                        AND Total_GWP IS NOT NULL
                        AND A1toA3_GWP != 0
                        """).fetchall()
emissions_reinf_values = [row[0] for row in emissions_reinf]

EPD_reinf = cursor.execute("""
                      SELECT Total_GWP FROM products
                      WHERE MATERIAL LIKE '%steel_reinforcing_bar%'
                      AND Total_GWP IS NOT NULL
                      AND "SOURCE" NOT LIKE '%KBOB%'
                      """).fetchall()
EPD_reinf_values = [row[0] for row in EPD_reinf]

KBOB_reinf = cursor.execute(
                            """
                            SELECT Total_GWP FROM products
                            WHERE MATERIAL LIKE '%steel_reinforcing_bar%'
                            AND "SOURCE" LIKE '%KBOB%'
                            """).fetchall()
KBOB_reinf_values = [row[0] for row in KBOB_reinf]

reinf_min = min(EPD_reinf_values)
reinf_max =max(EPD_reinf_values)

#------------------------------------------------------------------------------------------------------------------------
#extract values for steel
emissions_steel = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE MATERIAL LIKE '%structural_steel_profile%'
                        AND Total_GWP IS NOT NULL
                        AND Total_GWP != 0
                        """).fetchall()
emissions_steel_values = [row[0] for row in emissions_steel]

EPD_steel = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE MATERIAL LIKE '%structural_steel_profile%'
                        AND Total_GWP IS NOT NULL
                        AND "SOURCE" NOT LIKE '%KBOB%'
                        AND "SOURCE" NOT LIKE '%Stalia%'
                        """).fetchall()
EPD_steel_values = [row[0] for row in EPD_steel]

KBOB_steel = cursor.execute("""
                        SELECT A1toA3_GWP FROM products
                        WHERE MATERIAL LIKE '%structural_steel_profile%'
                        AND "SOURCE" LIKE '%KBOB%'
                        """).fetchall()
KBOB_steel_values = [row[0] for row in KBOB_steel]

#------------------------------------------------------------------------------------------------------------------------
#plot concrete
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 14,
    'legend.fontsize': 12,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,

    'font.family': 'Times New Roman',

    # Figure aesthetics
    'figure.facecolor': 'white',
    #'figure.autolayout': True,
    'legend.frameon': False,

    # Bar styling
    'patch.edgecolor': 'black',  # outline color for bars
    'patch.linewidth': 0.8,  # outline thickness
    'patch.facecolor': '#66c2a5',  # default bar color (optional)
    'patch.force_edgecolor': True  # ensures edgecolor is applied even if alpha is set

})

title = ['Beton', 'Holz', 'Bewehrungsstahl', 'Baustahl']
ymin = 0
ymax = 13
xmin = 0
xmax = 1400


import seaborn as sns
color = sns.color_palette()

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))#, sharey=True)

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
ax1.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color=color[2], alpha=0.4, label='EPD Histogramm')
ax1.scatter(EPD_concrete_values, np.random.normal(2, 0.5, len(EPD_concrete_values)), alpha=0.7, label='EPD Einzelwerte', color = 'white', edgecolor='black')

ax1.set_yticks(range(ymin, ymax))
ax1.set_ylim(ymin,ymax)
ax1.set_xlim(xmin,xmax)

ax1.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax1.set_ylabel('# EPD$_{Beton}$')
ax1.set_title(title[0])

# Add vertical lines and text for KBOB values
ax1.axvline(KBOB_concrete_values[0], linestyle='--', alpha=0.9, color=color[3], label='KBOB Hochbaubeton')

ax1.axvline(statistics.mean(EPD_concrete_values), linestyle='-', alpha=0.9, color=color[2], label= r'Mittelwert $\mu$')
ax1.fill_betweenx([ymin, ymax], np.quantile(EPD_concrete_values,0.1), np.quantile(EPD_concrete_values,0.9), color = color[7], alpha=0.1, zorder=0)
ax1.text(np.quantile(EPD_concrete_values, 0.1)-20, ymax-1, '80%', color = color[7], alpha=0.8)

ax1.legend(loc='upper right')

#------------------------------------------plot wood------------------------------------------------------------------
# Combine all datasets to determine the range for bins
all_data_timber = EPD_timber_values #+ Ecoinvent_timber_values
EPD_timber_values_pos = [x for x in EPD_timber_values if x >= 0]

# Determine the bins based on all data
bins = np.histogram_bin_edges(EPD_timber_values_pos, bins='auto')
# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histograms
hist_EPD, _ = np.histogram(EPD_timber_values_pos, bins=bins)
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
ax2.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color=color[1], alpha=0.4, label='EPD Histogramm')
ax2.scatter(EPD_timber_values_pos, np.random.normal(2, 0.5, len(EPD_timber_values_pos)), alpha=0.7, label='EPD Emissions', color = 'white', edgecolor='black')

ax2.set_yticks(range(ymin, ymax))
ax2.set_ylim(ymin,ymax)
ax2.set_xlim(xmin,xmax)

ax2.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax2.set_ylabel('# EPD$_{Holz}$')
ax2.set_title(title[1])

# Add vertical lines and text for KBOB values
ax2.axvline(253, linestyle='--', alpha= 1, color=sns.color_palette("Paired")[1], label='KBOB BSH CH')
ax2.axvline(335, linestyle='-.', alpha=1, color=sns.color_palette("Paired")[0], label='KBOB BSH')
ax2.axvline(281, linestyle='--', alpha=1, color=sns.color_palette("Paired")[5], label='KBOB BSP CH')
ax2.axvline(395, linestyle='-.', alpha=1, color=sns.color_palette("Paired")[4], label='KBOB BSP')
ax2.axvline(288, linestyle='-.', alpha=1, color=sns.color_palette("Paired")[3], label='KBOB KVH')

ax2.axvline(statistics.mean(EPD_timber_values_pos), linestyle='-', alpha=0.9, color=color[1], label = r'Mittelwert $\mu$')
ax2.fill_betweenx([ymin, ymax], np.quantile(EPD_timber_values_pos,0.1), np.quantile(EPD_timber_values_pos,0.9), color = color[7], alpha=0.1, zorder=0)
ax2.text(np.quantile(EPD_timber_values_pos, 0.1)-3, ymax-1, '80%', color = color[7], alpha=0.9)

ax2.legend(loc='upper right')

#----------------------------------------------plot Betonstahl----------------------------------------------------------
# Combine all datasets to determine the range for bins
all_data_reinf = EPD_reinf_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_reinf, bins='auto')

# Calculate the bin centers
bin_centers_EPD = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histogram for EPD_timber_values
hist_EPD, _ = np.histogram(EPD_reinf_values, bins=bins)

# Create the scatter plot for EPD values
bin_widths = np.diff(bins)
ax3.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color=color[0], alpha=0.4,label = 'EPD Histogramm')
ax3.scatter(EPD_reinf_values, np.random.normal(2, 0.5, len(EPD_reinf_values)), alpha=0.7, label='EPD Einzelwerte', color = 'white', edgecolor='black')

ax3.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax3.set_ylabel('# EPD$_{Bewehrung}$')
ax3.set_title('Betonstahl')

ax3.set_yticks(range(ymin, ymax))
ax3.set_ylim(ymin,ymax)
ax3.set_xlim(xmin,xmax)

# Add vertical lines and text for KBOB values
ax3.axvline(KBOB_reinf_values[0], linestyle='--', alpha = 0.7, color=color[3], label='KBOB Armierungsstahl')

ax3.axvline(368, color='grey',linestyle='dotted')
ax3.text(368+10, 4, 'Stahl Gerlafingen', fontsize = 12, rotation=90, color='grey')

# Add vertical lines and text for KBOB values
ax3.axvline(statistics.mean(EPD_reinf_values), linestyle='-', alpha=0.9, color=color[0], label=r'Mittelwert $\mu$')
ax3.fill_betweenx([ymin, ymax], np.quantile(EPD_reinf_values,0.1), np.quantile(EPD_reinf_values,0.9), color = color[7], alpha=0.1, zorder=0)
ax3.text(np.quantile(EPD_reinf_values, 0.1)+10, ymax-1, '80%', color = color[7], alpha=0.9)

ax3.legend(loc='upper right')

#------------------------------------------------------------------------------------------------------------------------
#plot steel
# Combine all datasets to determine the range for bins
all_data_steel = EPD_steel_values

# Determine the bins based on all data
bins = np.histogram_bin_edges(all_data_steel, bins=8)

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])
# Calculate the bin centers
bin_centers_EPD = 0.5 * (bins[:-1] + bins[1:])

# Calculate the histogram for EPD_timber_values
hist_EPD, _ = np.histogram(EPD_steel_values, bins=bins)

# Create the scatter plot for EPD values
bin_widths = np.diff(bins)
ax4.bar(bin_centers_EPD, hist_EPD, width=bin_widths[:len(bin_centers_EPD)], color=color[9], alpha=0.3, label = 'EPD Histogramm')
ax4.scatter(EPD_steel_values, np.random.normal(2, 0.5, len(EPD_steel_values)), alpha=0.7, label='EPD Einzelwerte', color = 'white', edgecolor='black')

ax4.set_yticks(range(ymin, ymax))
ax4.set_ylim(ymin,ymax)
ax4.set_xlim(xmin,xmax)

ax4.set_xlabel('Total GWP [kg CO$_2$-eq/t]')
ax4.set_ylabel('# EPD$_{Baustahl}$')
ax4.set_title('Baustahl')

# Add vertical lines and text for KBOB values
ax4.axvline(KBOB_steel_values[0], linestyle='--', alpha = 0.7, color=color[3], label = 'KBOB Baustahl')

# Add vertical lines and text for KBOB values
ax4.axvline(statistics.mean(EPD_steel_values), linestyle='-', alpha=0.9, color=color[9], label=r'Mittelwert $\mu$')
ax4.fill_betweenx([ymin, ymax], np.quantile(EPD_steel_values,0.1), np.quantile(EPD_steel_values,0.9), color = color[7], alpha=0.1, zorder=0)
ax4.text(np.quantile(EPD_steel_values, 0.1)+10, ymax-1, '80%', color = color[7], alpha=0.9)

ax4.legend(loc='upper right')

plt.show()

