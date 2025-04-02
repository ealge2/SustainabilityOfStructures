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
section = struct_analysis.RibbedConcrete(concrete1, rebar1, 16, 5.0, 0.40, 1.20, 0.18, 0.012, 0.15, 0.012, 0.15, 0.03, 7, 0.01, 0.15, 2, 2, 0.03 )

# create floor structure for solid wooden cross-section
bodenaufbau = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                                 ["'Unterlagsboden Zement, 85 mm'", False, False], ["'Glaswolle'", 0.03, False]]
print(type(bodenaufbau))

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
print("Querschnittswerte:")
print("beff = ", round(section.b_eff,5), "m")
print("zs = ", round(section.z_s,5), "[m]")
print("Iy = ", round(section.iy,5), "[m4]")
print("d_PB = ", section.d_PB, "[m]")
print("d_slab = ", section.d, "[m]")
print ("fsd= ", section.rebar_type.fsd)
print ("fcd = ", section.concrete_type.fcd)
print ("fctm = ", section.concrete_type.fctm)
print("As = ", section.as_PB_p)

print("Mu_max_slab = ", round(section.mu_max_slab), "[Nm]")
print("Vu_slab_p = ", round(section.vu_p), "[N]")
print("Mr = ",  round(section.mr_p), "[Nm]")
print("Mu_PB_max = ", round(section.mu_max), "[Nm]")
print("Mu_PB_min = ", round(section.mu_min), "[Nm]")
print("Mr_PB = ",  round(section.mr_pb_p), "[Nm]")
print("Vu_PB_p = ", round(section.vu_PB_p), "[N]")
print("Vu_PB_n = ", round(section.vu_PB_n), "[N]")

print("g0k= ", section.g0k)

print()
print("Querschnittsklasse", section.qs_class_p, section.qs_class_n, section.qs_class_p_slab, section.qs_class_n_slab)
print("qu = ", member.qu, "[N/m]")
print("Vu+ = ", round(member.section.vu_p), "[N], Vu- = ", round(member.section.vu_n), "[N]")

member.calc_qk_zul_gzt()
print("qk_zul_GZT = ", round(member.qk_zul_gzt), "[N/m]")
print("w,inst,adm = ", member.w_install_adm)
print("lw_install= ", member.requirements.lw_install)
print("w,use,adm = ", member.w_use_adm)
print("lw_use= ", member.requirements.lw_use)
print("w,app,adm = ", member.w_app_adm)
print("lw_app= ", member.requirements.lw_app)

print("Feuerwiderstand:")
member.get_fire_resistance()
print(member.fire_resistance)

print("qrare =", member.q_rare)
print("qpers =", member.q_per)
print("qfreq =", member.q_freq)

print("w_install = ", round(member.w_install,5))
print("w_install_ger = ", round(member.w_install_ger,5))
print("w_use = ", round(member.w_use,5))
print("w_use_ger = ", round(member.w_use_ger,5))
print("w_app = ", round(member.w_app,5))
print("w_app_ger = ", round(member.w_app_ger,5))




print("1. EF = ", member.f1)

