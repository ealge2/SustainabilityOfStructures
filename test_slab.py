# file contains code for verification of members with wooden cross-sections


# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization_2D
import struct_optimization  # file with code for structural optimization
#import matplotlib.pyplot as plt

# INPUT
# create dummy-database
database_name = "database_250702.db"  # define database name
#create_dummy_database.create_database(database_name)  # create database

# create material for reinforced concrete cross-section, derive corresponding design values
concrete1 = struct_analysis.ReadyMixedConcrete("'C25/30'", database_name)  # create a Wood material object
concrete1.get_design_values()
rebar1 = struct_analysis.SteelReinforcingBar("'B500B'", database_name)  # create a Wood material object
rebar1.get_design_values()

# create reinforced concrete rectangular cross-section
section = struct_analysis.RectangularConcrete(concrete1, rebar1, 1.0, 0.24, 0.012, 0.15, 0.012, 0.15, 0.0, 0.15, 0)

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

length_x = 8
length_y = 8
support = "LL-frei"

# create slab system
system = struct_analysis.Slab(length_x,length_y,support)


# create rc member
member = struct_analysis.Member2D(section, system, bodenaufbau_rc, requirements, g2k, qk)
opt_section = struct_optimization_2D.get_optimized_section(member, "ENV", "GWP", 25)
print("opt section = ", opt_section.h)

print("d =", section.d)
print("mu_max= ", round(section.mu_max,2))
print("alpha_m: ",system.alpha_m)
print("alpha_m_neg: ",system.alpha_m_neg)
print("mu_min= ", round(section.mu_min,2))


print()
print("as = ", opt_section.bw)


print("mr_p =", section.mr_p)
print("x/d =", section.x_p/section.d)
print()
print("qu =", round(member.qu,2))
print("vu = ", member.section.vu_p, member.section.vu_n)

member.calc_qk_zul_gzt()
print("qk_zul_gzt =", member.qk_zul_gzt)
print("Feuerwiderstand:")
member.get_fire_resistance()
print(member.fire_resistance)


print("w_inst_adm=", round(member.w_install_adm,5))
print("w_use_adm=", round(member.w_use_adm,5))
print("w_app_adm=", round(member.w_app_adm,5))