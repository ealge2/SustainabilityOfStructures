# file contains code for generating fundamental plots of first example (simple supported beam)

# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import matplotlib.pyplot as plt
import numpy as np

# INPUT
# create dummy-database
# database_name = "dummy_sustainability.db"  # define database name
# create_dummy_database.create_database(database_name)  # create database

# define database
database_name = "database_250514.db"

# create material for wooden cross-section, derive corresponding design values
timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber1.get_design_values()


# create materials for reinforced concrete cross-section, derive corresponding design values
concrete1 = struct_analysis.ReadyMixedConcrete("'C25/30'", database_name)
concrete1.get_design_values()
reinfsteel1 = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
reinfsteel1.get_design_values()



# Placeholder values for demonstration
materials = ['Wood', 'Concrete', 'Reinforcing Steel']
Strength = [timber1.fmd, concrete1.fcd, reinfsteel1.fsd]  # in MPa
Strength_array = np.array(Strength)
E_modulus = [timber1.Emmean, concrete1.Ecm, reinfsteel1.Es]  # in MPa
E_modulus_array = np.array(E_modulus)
GWP = [timber1.GWP, concrete1.GWP, reinfsteel1.GWP]  # in kg CO2-eq/kg
GWP_array = np.array(GWP)
ratioStiffness = GWP_array/E_modulus_array
ratioStrength = GWP_array/Strength_array

# Create figure and axis
fig, ax1 = plt.subplots()

# Bar plot for E-modulus
color = 'tab:blue'
ax1.set_xlabel('Material')
ax1.set_ylabel('GWP/ExA (kgCO2-eq/N*m)', color=color)
bars = ax1.bar(materials, ratioStiffness, color=color, alpha=0.6, label='E-Modulus')
ax1.tick_params(axis='y', labelcolor=color)

"""
# Create a second y-axis for GWP
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('GWP (kg CO2-eq/kg)', color=color)
line = ax2.plot(materials, GWP, color=color, marker='o', label='GWP')
ax2.tick_params(axis='y', labelcolor=color)
"""
# Title and layout
plt.title('Comparison of E-Modulus and GWP for Different Materials')
fig.tight_layout()

# Show plot
plt.show()


# Create figure and axis
fig, ax1 = plt.subplots()

# Bar plot for E-modulus
color = 'tab:blue'
ax1.set_xlabel('Material')
ax1.set_ylabel('GWP/fyA (kgCO2-eq/N*m)', color=color)
bars = ax1.bar(materials, ratioStrength, color=color, alpha=0.6, label='E-Modulus')
ax1.tick_params(axis='y', labelcolor=color)

"""
# Create a second y-axis for GWP
ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('GWP (kg CO2-eq/kg)', color=color)
line = ax2.plot(materials, GWP, color=color, marker='o', label='GWP')
ax2.tick_params(axis='y', labelcolor=color)
"""
# Title and layout
plt.title('Comparison of Strength and GWP for Different Materials')
fig.tight_layout()

# Show plot
plt.show()

