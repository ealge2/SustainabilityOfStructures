# file contains code for verification of members with wooden cross-sections


# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import matplotlib.pyplot as plt

# INPUT
# create dummy-database
database_name = "dummy_sustainability.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# create material for wooden cross-section, derive corresponding design values
timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber1.get_design_values()

# create initial wooden rectangular cross-section
section = struct_analysis.RectangularWood(timber1, 1.0, 0.24)

# create floor structure for solid wooden cross-section
bodenaufbau = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                                 ["'Unterlagsboden Zement, 85 mm'", False, False], ["'Glaswolle'", 0.03, False],
                                 ["'Kies gebrochen'", 0.06, False]]
bodenaufbau_wd = struct_analysis.FloorStruc(bodenaufbau, database_name)

requirements = struct_analysis.Requirements()

# define loads on member
g2k = 0.75e3  # n.t. Einbauten
qk = 2e3  # Nutzlast

# define service limit state criteria
req = struct_analysis.Requirements()

# define system length
length = 5.8

# create simple supported beam system
system = struct_analysis.BeamSimpleSup(length)

# create wooden member
member = struct_analysis.Member1D(section, system, bodenaufbau_wd, requirements, g2k, qk)

print(member.f1)
print(member.a_ed)  #XXXXXXXXXXXXXXXXX result of a_ed looks wrong. has to be corrected.XXXXXXXXXXXXXXXX
