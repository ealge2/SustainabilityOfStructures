# file contains code for verification of members with wooden cross-sections


# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
#import struct_optimization  # file with code for structural optimization
#import matplotlib.pyplot as plt

# INPUT
# create dummy-database
database_name = "dummy_sustainability_1.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# create material for reinforced concrete cross-section, derive corresponding design values
concrete1 = struct_analysis.ReadyMixedConcrete("'C25/30'", database_name)  # create a Wood material object
concrete1.get_design_values()
rebar1 = struct_analysis.SteelReinforcingBar("'B500B'", database_name)  # create a Wood material object
rebar1.get_design_values()

# create reinforced concrete rectangular cross-section
section = struct_analysis.RibbedConcrete(concrete1, rebar1, 16, 5.0, 0.40, 1.20, 0.18, 0.012, 0.15, 0.012, 0.15, 0.03, 7, 0, 0.15, 0.01, 0.15, 2 )

# create floor structure for solid wooden cross-section
bodenaufbau = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                                 ["'Unterlagsboden Zement, 85 mm'", False, False], ["'Glaswolle'", 0.03, False]]
bodenaufbau_rc = struct_analysis.FloorStruc(bodenaufbau, database_name)

requirements = struct_analysis.Requirements()

# define loads on member
g2k = 0.75e3  # n.t. Einbauten
qk = 2e3  # Nutzlast

# define service limit state criteria
req = struct_analysis.Requirements()

# define system length
length = 16

# create simple supported beam system
system = struct_analysis.BeamSimpleSup(length)

# create wooden member
member = struct_analysis.Member1D(section, system, bodenaufbau_rc, requirements, g2k, qk)
print("beff = ", section.b_eff, "m")
print("Iy = ", section.iy, "[m4]")
print("zs = ", section.z_s, "[m]")
print(" d = ", section.d_PB, "[m]")
print ("fsd= ", section.rebar_type.fsd)
print ("fcd = ", section.concrete_type.fcd)
print("As = ", section.as_PB_p)
print("Mu_max = ", round(section.mu_max), "[Nm]")
print("Mr = ",  round(section.mr_p), "[Nm]")
print("Mu_PB_max = ", round(section.mu_PB_max), "[Nm]")
print("Mr_PB = ",  round(section.mr_pb_p), "[Nm]")


print("x/d = ", round(section.x_p/section.d*1000,1), "[mm]")
print()
print("qu = ", member.qu, "[N/m]")
print("Vu+ = ", round(member.section.vu_p), "[N], Vu- = ", round(member.section.vu_n), "[N]")

member.calc_qk_zul_gzt()
print("qk_zul_GZT = ", round(member.qk_zul_gzt), "[N/m]")
print("Feuerwiderstand:")
member.get_fire_resistance()
print(member.fire_resistance)
