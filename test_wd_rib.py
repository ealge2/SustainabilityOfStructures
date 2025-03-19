# file contains code for verification of members with wooden cross-sections


# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
#import matplotlib.pyplot as plt

# INPUT
# create dummy-database
database_name = "dummy_sustainability_1.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# create material for wooden cross-section, derive corresponding design values
timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber1.get_design_values()

# create wooden rectangular cross-section
section = struct_analysis.RibWood(timber1, timber1, timber1, 6.25,0.12, 0.26, 0.625, 0.027, 0.027)


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

print("bc_ef = ", section.bc_ef, "[m]")
print("bt_ef = ", section.bt_ef, "[m]")
print("n =", section.n)
print("n_inf =", section.n_inf)
print("z_s = ", section.z_s)
print("iy= ", section.iy)
print("iy_inf= ", section.iy_inf)

print("g1k  =", member.g1k)
print("g0k =", member.g0k)
print("gk =", member.gk)

print("ei =", round(section.ei1/.625,5))
print("E =", section.wood_type_1.Emmean)
print("E =", section.wood_type_2.Emmean)
print("E =", section.wood_type_3.Emmean)

print("f1 =", member.f1, "[Hz]")
print(member.a_ed)
print(section.mu_max)
print(member.section.vu_p, member.section.vu_n)
member.calc_qk_zul_gzt()
print(member.qk_zul_gzt)
print("Feuerwiderstand:")
member.get_fire_resistance()
print(member.fire_resistance)
print(" ")