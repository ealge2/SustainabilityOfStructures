import matplotlib.pyplot as plt
import sqlite3
import numpy as np
import statistics
import struct_analysis  # file with code for structural analysis
import os

# Define the relative path to the database file
database_path = os.path.join("Database", "database_250326.db")
# Use this path to connect to the database
connection = sqlite3.connect(database_path)
# create cursor object
cursor = connection.cursor()
"""
# create or open database sustainability
connection = sqlite3.connect('database_250326.db')
# create cursor object
cursor = connection.cursor()
"""
#------------------------------------------------------------------------------------------------------------------------
#extract values for concrete
#

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
Density_concrete = cursor.execute(
                    """
                    SELECT DENSITY FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%Betonsortenrechner%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND "source [string]" NOT LIKE '%KBOB%'
                    """).fetchall()
Density_EPDconcrete_values = [row[0] for row in Density_concrete]

EPD_concreteC2530 = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND MECH_PROP LIKE '%C25/30%'
                    AND "source [string]" NOT LIKE '%Betonsortenrechner%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND "source [string]" NOT LIKE '%KBOB%'
                    """).fetchall()
EPD_concreteC2530_values = [row[0] for row in EPD_concreteC2530]

Density_concreteC2530 = cursor.execute(
                    """
                    SELECT DENSITY FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND MECH_PROP LIKE '%C25/30%'
                    AND "source [string]" NOT LIKE '%Betonsortenrechner%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND "source [string]" NOT LIKE '%KBOB%'
                    """).fetchall()
Density_EPDconcreteC2530_values = [row[0] for row in Density_concreteC2530]

EPD_concreteC3037 = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND MECH_PROP LIKE '%C30/37%'
                    AND "source [string]" NOT LIKE '%Betonsortenrechner%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND "source [string]" NOT LIKE '%KBOB%'
                    """).fetchall()
EPD_concreteC3037_values = [row[0] for row in EPD_concreteC3037]

Density_concreteC3037 = cursor.execute(
                    """
                    SELECT DENSITY FROM products
                    WHERE type LIKE '%ready mixed concrete%'
                    AND A1toA3_GWP IS NOT NULL
                    AND MECH_PROP LIKE '%C30/37%'
                    AND "source [string]" NOT LIKE '%Betonsortenrechner%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND "source [string]" NOT LIKE '%KBOB%'
                    """).fetchall()
Density_EPDconcreteC3037_values = [row[0] for row in Density_concreteC3037]

#------------------------------------------------------------------------------------------------------------------------
#extract values for wood

EPD_timber_GL24h = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%KBOB%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND MECH_PROP LIKE '%GL24h%'
                    AND EPD_ID NOT LIKE '%verifizierung%'
                    """).fetchall()
EPD_timber_GL24h_values = [row[0] for row in EPD_timber_GL24h]

Density_timber_GL24h = cursor.execute(
                    """
                    SELECT DENSITY FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%KBOB%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND MECH_PROP LIKE '%GL24h%'
                    AND EPD_ID NOT LIKE '%verifizierung%'
                    """).fetchall()
Density_timber_GL24h_values = [row[0] for row in Density_timber_GL24h]

EPD_timber_C24 = cursor.execute(
                    """
                    SELECT Total_GWP FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%KBOB%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND MECH_PROP LIKE '%C24%'
                    AND EPD_ID NOT LIKE '%verifizierung%'
                    """).fetchall()
EPD_timber_C24_values = [row[0] for row in EPD_timber_C24]

Density_timber_C24 = cursor.execute(
                    """
                    SELECT DENSITY FROM products
                    WHERE material_0 LIKE '%timber%'
                    AND Total_GWP IS NOT NULL
                    AND "source [string]" NOT LIKE '%KBOB%'
                    AND "source [string]" NOT LIKE '%Ecoinvent%'
                    AND MECH_PROP LIKE '%C24%'
                    AND EPD_ID NOT LIKE '%verifizierung%'
                    """).fetchall()
Density_timber_C24_values = [row[0] for row in Density_timber_C24]


#------------------------------------------------------------------------------------------------------------------------
#extract values for reinforcement

EPD_reinf = cursor.execute("""
                      SELECT Total_GWP FROM products
                      WHERE type LIKE '%steel reinforcing bar%'
                      AND Total_GWP IS NOT NULL
                      AND "source [string]" NOT LIKE '%KBOB%'
                      """).fetchall()
EPD_reinf_values = [row[0] for row in EPD_reinf]

Density_reinf = cursor.execute("""
                      SELECT Density FROM products
                      WHERE type LIKE '%steel reinforcing bar%'
                      AND Total_GWP IS NOT NULL
                      AND "source [string]" NOT LIKE '%KBOB%'
                      """).fetchall()
Density_reinf_values = [row[0] for row in Density_reinf]

reinf_min = min(EPD_reinf_values)
reinf_max =max(EPD_reinf_values)

#------------------------------------------------------------------------------------------------------------------------
#extract values for steel

EPD_steel = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE type LIKE '%structural steel profile%'
                        AND Total_GWP IS NOT NULL
                        AND "source [string]" NOT LIKE '%KBOB%'
                        """).fetchall()
EPD_steel_values = [row[0] for row in EPD_steel]

#------------------------------------------------------------------------------------------------------------------------
#extract values for prestressingSteel

EPD_prestressingSteel = cursor.execute("""
                        SELECT Total_GWP FROM products
                        WHERE type LIKE '%prestressing steel%'
                        AND Total_GWP IS NOT NULL
                        AND "source [string]" NOT LIKE '%KBOB%'
                        """).fetchall()
EPD_prestressingSteel_values = [row[0] for row in EPD_prestressingSteel]

#------------------------------------------------------------------------------------------------------------------------
# define database
database_name = "database_250514.db"
# create materials for reinforced concrete cross-section, derive corresponding design values
concrete1 = struct_analysis.ReadyMixedConcrete("'C30/37'", database_name)
concrete1.get_design_values()
reinfsteel1 = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
reinfsteel1.get_design_values()
# create material for wooden cross-section, derive corresponding design values
timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber1.get_design_values()

# Extract Mean Data
materials = ['Wood', 'Concrete', 'Reinforcing Steel']
timber_GL24h_fc0d = 14500000
timber_C24_fc0d = 12400000
timber_GL24h_ft0d = 12800000
timber_C24_ft0d = 8500000
timber_GL24h_Emmean = 11500000000
timber_C24_Emmean =  11000000000
concrete_C2530_fcd = 16500000
concrete_C3037_fcd = 20000000
concrete_C8590_fcd = 56700000
steel_fsd = 550000000/1.05
steel_Es = 210000000000
presteel_fsd = 0.85*1860000000/1.15
presteel_Es = 200000000000

Strength = [timber1.fmd, concrete1.fcd, reinfsteel1.fsd]  # in MPa
Strength_array = np.array(Strength)
E_modulus = [timber1.Emmean, concrete1.Ecm, reinfsteel1.Es]  # in MPa
E_modulus_array = np.array(E_modulus)
GWP = [timber1.GWP, concrete1.GWP, reinfsteel1.GWP]  # in kg CO2-eq/kg
GWP_array = np.array(GWP)

ratioStiffness = GWP_array/E_modulus_array
ratioStrength = GWP_array/Strength_array

#------------------------------------------------------------------------------------------------------------------------
# filter out negative values

EPD_concrete_values = np.array(EPD_concrete_values, dtype=float)
Density_EPDconcrete_values = np.array(Density_EPDconcrete_values, dtype=float)
non_zero_EPD_concrete = EPD_concrete_values > 0
EPD_concrete_values = EPD_concrete_values*Density_EPDconcrete_values
EPD_concrete_values = EPD_concrete_values[non_zero_EPD_concrete]

EPD_concreteC2530_values = np.array(EPD_concreteC2530_values, dtype=float)
Density_EPDconcreteC2530_values = np.array(Density_EPDconcreteC2530_values, dtype=float)
non_zero_EPD_concreteC2530 = EPD_concreteC2530_values > 0
EPD_concreteC2530_values = EPD_concreteC2530_values*Density_EPDconcreteC2530_values
EPD_concreteC2530_values = EPD_concreteC2530_values[non_zero_EPD_concreteC2530]

EPD_concreteC3037_values = np.array(EPD_concreteC3037_values, dtype=float)
Density_EPDconcreteC3037_values = np.array(Density_EPDconcreteC3037_values, dtype=float)
non_zero_EPD_concreteC3037 = EPD_concreteC3037_values > 0
EPD_concreteC3037_values = EPD_concreteC3037_values*Density_EPDconcreteC3037_values
EPD_concreteC3037_values = EPD_concreteC3037_values[non_zero_EPD_concreteC3037]

EPD_timber_GL24h_values = np.array(EPD_timber_GL24h_values, dtype=float)
Density_timber_GL24h_values = np.array(Density_timber_GL24h_values, dtype=float)
non_zero_EPD_GL24h_timber = EPD_timber_GL24h_values > 0
EPD_timber_GL24h_values = EPD_timber_GL24h_values*Density_timber_GL24h_values
EPD_timber_GL24h_values = EPD_timber_GL24h_values[non_zero_EPD_GL24h_timber]

EPD_timber_C24_values = np.array(EPD_timber_C24_values, dtype=float)
Density_timber_C24_values = np.array(Density_timber_C24_values, dtype=float)
non_zero_EPD_C24_timber = EPD_timber_C24_values > 0
EPD_timber_C24_values = EPD_timber_C24_values*Density_timber_C24_values
EPD_timber_C24_values = EPD_timber_C24_values[non_zero_EPD_C24_timber]

EPD_reinf_values = np.array(EPD_reinf_values)
non_zero_EPD_reinf = EPD_reinf_values > 0
EPD_reinf_values = EPD_reinf_values*7850
EPD_reinf_values = EPD_reinf_values[non_zero_EPD_reinf]

EPD_steel_values = np.array(EPD_steel_values)
non_zero_EPD_steel = EPD_steel_values > 0
EPD_steel_values = EPD_steel_values*7850
EPD_steel_values = EPD_steel_values[non_zero_EPD_steel]

EPD_prestressingSteel_values = np.array(EPD_prestressingSteel_values)
non_zero_EPD_prestressingSteel = EPD_prestressingSteel_values > 0
EPD_prestressingSteel_values = EPD_prestressingSteel_values*7850
EPD_prestressingSteel_values = EPD_prestressingSteel_values[non_zero_EPD_prestressingSteel]

#------------------------------------------------------------------------------------------------------------------------
# PLOT DATA

# Normalize values by Stiffness
normalized_by_Stiffness = [
    np.array(EPD_concreteC3037_values) / concrete1.Ecm*1000000,
    np.array(EPD_timber_GL24h_values) / timber_GL24h_Emmean*1000000,
    np.array(EPD_timber_C24_values) / timber_C24_Emmean*1000000,
    np.array(EPD_reinf_values) / reinfsteel1.Es*1000000,
    np.array(EPD_steel_values) / steel_Es*1000000,
    np.array(EPD_prestressingSteel_values) / presteel_Es*1000000,
]


# Normalize values by respective compressive strength values
normalized_by_cstrength = [
    np.array(EPD_concreteC3037_values) / concrete_C3037_fcd*1000,
    np.array(EPD_timber_GL24h_values) / timber_GL24h_fc0d*1000,
    np.array(EPD_timber_C24_values) / timber_C24_fc0d*1000,
    np.array(EPD_reinf_values) / reinfsteel1.fsd*1000,
    np.array(EPD_steel_values) / steel_fsd*1000,
    np.array(EPD_prestressingSteel_values) / presteel_fsd*1000
]

# Normalize values by respective tensile strength values
normalized_by_tstrength = [
    np.array(EPD_timber_GL24h_values) / timber_GL24h_ft0d*1000,
    np.array(EPD_timber_C24_values) / timber_C24_ft0d*1000,
    np.array(EPD_reinf_values) / reinfsteel1.fsd*1000,
    np.array(EPD_steel_values) / steel_fsd*1000,
    np.array(EPD_prestressingSteel_values) / presteel_fsd*1000
]

# Volumetric
volumetricGWP = [
    np.array(EPD_concreteC3037_values),
    np.array(EPD_timber_GL24h_values),
    np.array(EPD_timber_C24_values),
    np.array(EPD_reinf_values),
    np.array(EPD_steel_values),
    np.array(EPD_prestressingSteel_values)
]
# ______________________________________________

plt.rcParams['font.family'] = 'Times New Roman'  # globale Schrift
plt.rcParams['font.size'] = 10                  # optionale Grundgröße
plt.rcParams['axes.titlesize'] = 10             # Titelgröße
plt.rcParams['axes.labelsize'] = 10             # Achsentitelgröße
plt.rcParams['xtick.labelsize'] = 10            # Tick-Größen
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10


# Create subplots english
fig, axs = plt.subplots(2, 2, figsize=(8, 5))
ax1, ax2, ax3, ax4 = axs.flatten()
# Left subplot: normalized by Stiffness
ax1.boxplot(normalized_by_Stiffness, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax1.set_ylim(bottom=0)   # untere Grenze auf 0
ax1.set_title('GWP normalised by stiffness')
ax1.set_xticklabels(['C37/30', 'GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax1.set_ylabel('GWP·\u03C1 / E [kgCO$_2$eq/(MN·m)]')

# Right subplot: normalized by strength
ax2.boxplot(normalized_by_cstrength, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax2.set_ylim(bottom=0)   # untere Grenze auf 0
ax2.set_title('GWP normalised by compressive strength')
ax2.set_xticklabels(['C37/30', 'GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax2.set_ylabel('GWP·\u03C1 / f [kgCO$_2$eq/(kN·m)]')

# Right subplot: normalized by strength
ax3.boxplot(normalized_by_tstrength, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax3.set_ylim(bottom=0)   # untere Grenze auf 0
ax3.set_title('GWP normalised by tensile strength')
ax3.set_xticklabels(['GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax3.set_ylabel('GWP·\u03C1 / f [kgCO$_2$eq/(kN·m)]')

# Right subplot: volumetric
ax4.boxplot(volumetricGWP, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax4.set_ylim(bottom=0)   # untere Grenze auf 0
ax4.set_title('volumetric GWP')
ax4.set_xticklabels(['C37/30', 'GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax4.set_ylabel('GWP·\u03C1 [kgCO$_2$eq/m$^3$]')
plt.tight_layout()
plt.show()

fig, ax1 = plt.subplots(1, 1, figsize=(4, 3))
ax1.boxplot(volumetricGWP[:3], patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax1.set_ylim(bottom=0)   # untere Grenze auf 0
ax1.set_title('volumetric GWP')
ax1.set_xticklabels(['C37/30', 'GL24h', 'C24'], rotation=45)
ax1.set_ylabel('GWP·\u03C1 [kgCO$_2$eq/m$^3$]')
plt.tight_layout()
plt.show()

# Deutsch
# Create subplots
fig, axs = plt.subplots(2, 2, figsize=(8, 5))
ax1, ax2, ax3, ax4 = axs.flatten()
# Left subplot: normalized by Stiffness
ax1.boxplot(normalized_by_Stiffness, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax1.set_ylim(bottom=0)   # untere Grenze auf 0
ax1.set_title('GWP bezogen auf die Steifigkeit')
ax1.set_xticklabels(['C37/30', 'GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax1.set_ylabel('GWP·\u03C1 / E [kgCO$_2$eq/(MN·m)]')

# Right subplot: normalized by strength
ax2.boxplot(normalized_by_cstrength, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax2.set_ylim(bottom=0)   # untere Grenze auf 0
ax2.set_title('GWP bezogen auf die Druckfestigkeit')
ax2.set_xticklabels(['C37/30', 'GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax2.set_ylabel('GWP·\u03C1 / f [kgCO$_2$eq/(kN·m)]')

# Right subplot: normalized by strength
ax3.boxplot(normalized_by_tstrength, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax3.set_ylim(bottom=0)   # untere Grenze auf 0
ax3.set_title('GWP bezogen auf die Zugfestigkeit')
ax3.set_xticklabels(['GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax3.set_ylabel('GWP·\u03C1 / f [kgCO$_2$eq/(kN·m)]')

# Right subplot: volumetric
ax4.boxplot(volumetricGWP, patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax4.set_ylim(bottom=0)   # untere Grenze auf 0
ax4.set_xticklabels(['C37/30', 'GL24h', 'C24','B500B', 'S 460', 'Y1860'], rotation=45)
ax4.set_ylabel('GWP·\u03C1 [kgCO$_2$eq/m$^3$]')
plt.tight_layout()
plt.show()

fig, ax1 = plt.subplots(1, 1, figsize=(4, 3))
ax1.boxplot(volumetricGWP[:3], patch_artist=True,
    boxprops=dict(facecolor='gray', color='gray'),  # Füllfarbe + Randfarbe
    medianprops=dict(color='black'),                # optional: Medianfarbe
    whiskerprops=dict(color='gray'),
    capprops=dict(color='gray'))
ax1.set_ylim(bottom=0)   # untere Grenze auf 0
ax1.set_xticklabels(['C37/30', 'GL24h', 'C24'], rotation=45)
ax1.set_ylabel('GWP·\u03C1 [kgCO$_2$eq/m$^3$]')
plt.tight_layout()
plt.show()
# Create subplots

fig, axs = plt.subplots(2, 2, figsize=(10, 10))
ax1, ax2, ax3, ax4 = axs.flatten()

# Left subplot: normalized by Stiffness
ax1.boxplot(normalized_by_Stiffness, patch_artist=True)
ax1.set_title('Normalized GWP by Stiffness')
ax1.set_xticklabels(['Concrete Average','C25/20', 'C85/90', 'Gl24h', 'C24','B500B'], rotation=45)
ax1.set_ylabel('GWP·\u03C1 / E [kgCO2eq/(MN·m)]')

# Right subplot: normalized by strength
ax2.boxplot(normalized_by_strength, patch_artist=True)
ax2.set_title('Normalized GWP by Compressive Strength')
ax2.set_xticklabels(['Concrete Average','C25/20', 'C85/90', 'Gl24h', 'C24','B500B'], rotation=45)
ax2.set_ylabel('GWP·\u03C1 / f [kgCO2eq/(kN·m)]')

# Right subplot: volumetric
ax3.boxplot(volumetricGWP, patch_artist=True)
ax3.set_title('Volumetric GWP')
ax3.set_xticklabels(['Concrete Average','C25/20', 'C85/90', 'Gl24h', 'C24','B500B'], rotation=45)
ax3.set_ylabel('GWP·\u03C1 [kgCO2eq/m3]')

# Right subplot: volumetric
ax4.boxplot(volumetricGWP[:5], patch_artist=True)
ax4.set_title('Volumetric GWP')
ax4.set_xticklabels(['Concrete Average','C25/20', 'C85/90', 'Gl24h', 'C24'], rotation=45)
ax4.set_ylabel('GWP·\u03C1 [kgCO2eq/m3]')

plt.tight_layout()
plt.show()

# Compute normalized values
"""
Ecoinvent_concrete_values = np.random.normal(300, 50, 30)
Betonsortenrechner_concrete_values = np.random.normal(280, 40, 30)
EPD_concrete_values = np.random.normal(290, 45, 30)
EPD_timber_values = np.random.normal(150, 30, 30)
EPD_reinf_values = np.random.normal(500, 60, 30)
"""

y_epd_conc_f = np.array(concrete1.fcd) / np.array(EPD_concrete_values) /1000
y_epd_conc_E = np.array(concrete1.Ecm) / np.array(EPD_concrete_values)/1000000
x_epd_conc = np.array(EPD_concrete_values)/1000

y_epd_concC2530_f = np.array(concrete_C2530_fcd) / np.array(EPD_concreteC2530_values)/1000
y_epd_concC2530_E = np.array(concrete1.Ecm) / np.array(EPD_concreteC2530_values)/1000000
x_epd_concC2530 = np.array(EPD_concreteC2530_values)/1000

y_epd_concC3037_f = np.array(concrete_C3037_fcd) / np.array(EPD_concreteC3037_values)/1000
y_epd_concC3037_E = np.array(concrete1.Ecm) / np.array(EPD_concreteC3037_values)/1000000
x_epd_concC3037 = np.array(EPD_concreteC3037_values)/1000

y_epd_timberGl24h_f = np.array(timber_GL24h_fc0d) / np.array(EPD_timber_GL24h_values)/1000
y_epd_timberGl24h_E = np.array(timber_GL24h_Emmean) / np.array(EPD_timber_GL24h_values)/1000000
x_epd_timberGl24h = np.array(EPD_timber_GL24h_values)/1000

y_epd_timberC24_f = np.array(timber_C24_fc0d) / np.array(EPD_timber_C24_values)/1000
y_epd_timberC24_E = np.array(timber_C24_Emmean) / np.array(EPD_timber_C24_values)/1000000
x_epd_timberC24 = np.array(EPD_timber_C24_values)/1000

y_epd_reinf_f = np.array(reinfsteel1.fsd) / np.array(EPD_reinf_values)/1000
y_epd_reinf_E = np.array(reinfsteel1.Es) / np.array(EPD_reinf_values)/1000000
x_epd_reinf = np.array(EPD_reinf_values)/1000


# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

# Left subplot: E/GWP vs GWP
ax1.scatter(x_epd_conc, y_epd_conc_E, label='Concrete EPD', alpha=0.7)
ax1.scatter(x_epd_concC2530, y_epd_concC2530_E, label='C25/30 EPD', alpha=0.7)
ax1.scatter(x_epd_concC3037, y_epd_concC3037_E, label='C85/90 EPD', alpha=0.7)
ax1.scatter(x_epd_timberGl24h, y_epd_timberGl24h_E, label='Gl24h', alpha=0.7)
ax1.scatter(x_epd_timberC24, y_epd_timberC24_E, label='C24', alpha=0.7)
ax1.scatter(x_epd_reinf, y_epd_reinf_E, label='Reinforcement B500B', alpha=0.7)
ax1.set_xlabel('GWP·\u03C1 [tCO2-eq/m3] ')
ax1.set_ylabel('E/ GWP·\u03C1 [MN·m / kg-CO2-eq]')
ax1.set_title('Stiffness')
ax1.legend()
ax1.grid(True)

# Right subplot: f/GWP vs GWP
ax2.scatter(x_epd_conc, y_epd_conc_f, label='Concrete EPD', alpha=0.7)
ax2.scatter(x_epd_concC2530, y_epd_concC2530_f, label='C25/30 EPD', alpha=0.7)
ax2.scatter(x_epd_concC3037, y_epd_concC3037_f, label='C85/90 EPD', alpha=0.7)
ax2.scatter(x_epd_timberGl24h, y_epd_timberGl24h_f, label='Gl24h', alpha=0.7)
ax2.scatter(x_epd_timberC24, y_epd_timberC24_f, label='C24', alpha=0.7)
ax2.scatter(x_epd_reinf, y_epd_reinf_f, label='Reinforcement B500B', alpha=0.7)
ax2.set_xlabel('GWP·\u03C1 [tCO2-eq/m3] ')
ax2.set_ylabel('f/ GWP·\u03C1 [kN·m / kg-CO2-eq]')
ax2.set_title('Compressive Strength')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()
