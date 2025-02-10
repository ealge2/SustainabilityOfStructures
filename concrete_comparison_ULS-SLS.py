# file contains code for generating fundamental plots of first example (simple supported beam)

# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import matplotlib.pyplot as plt

# INPUT
# create dummy-database
database_name = "dummy_sustainability.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# create materials for reinforced concrete cross-section, derive corresponding design values
concrete1 = struct_analysis.ReadyMixedConcrete("'C25/30'", database_name)
concrete1.get_design_values()
reinfsteel1 = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
reinfsteel1.get_design_values()

section_rc0 = struct_analysis.RectangularConcrete(concrete1, reinfsteel1, 1.0, 0.12, 0.014, 0.15, 0.01, 0.15)

bodenaufbau_rcdecke = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                       ["'Unterlagsboden Zement, 85 mm'", False, False],
                       ["'Glaswolle'", 0.03, False]]
bodenaufbau_rc = struct_analysis.FloorStruc(bodenaufbau_rcdecke, database_name)

# define loads on member
g2k = 0.75e3  # n.t. Einbauten
qk = 2e3  # Nutzlast

# define service limit state criteria
req = struct_analysis.Requirements()

# define system lengths for plot
lengths = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

#  define content of plot
to_plot = [[section_rc0, bodenaufbau_rc]]
criteria = ["ULS", "SLS1", "SLS2", "FIRE", "ENV"]
optima = ["GWP"]
plotted_data = [["h_struct", "[m]"], ["h_tot", "[m]"], ["GWP_struct", "[kg-CO2-eq]"], ["GWP_tot", "[kg-CO2-eq]"],
                ["cost_struct", "[CHF]"]]

# ANALYSIS
# max. number of iterations per optimization. Higher value leads to better results
max_iter = 100
member_list = []
legend = []
# create plot data
for i in to_plot:
    for criterion in criteria:
        for optimum in optima:
            members = []
            for length in lengths:
                sys = struct_analysis.BeamSimpleSup(length)
                member0 = struct_analysis.Member1D(i[0], sys, i[1], req, g2k, qk)
                opt_section = struct_optimization.get_optimized_section(member0, criterion, optimum, max_iter)
                opt_member = struct_analysis.Member1D(opt_section, sys, i[1], req, g2k, qk)
                members.append(opt_member)
            member_list.append(members)
            legend.append([i[0].section_type, criterion, optimum])

# plot figures
plt.figure(1)
data_max = [0, 0, 0, 0, 0, 0]
for i, members in enumerate(member_list):
    plotdata = [[], [], [], [], []]
    for mem in members:
        plotdata[0].append(mem.section.h)
        plotdata[1].append(mem.section.h + mem.floorstruc.h)
        plotdata[2].append(mem.section.co2)
        plotdata[3].append(mem.section.co2 + mem.floorstruc.co2)
        plotdata[4].append(mem.section.cost)
    sec_typ, cri, opt = legend[i]
    # set line color
    if sec_typ == "rc_rec":
        color = "tab:green"  # color for reinforced concrete
    elif sec_typ == "wd_rec":
        color = "tab:brown"  # color for wood
    else:
        color = "k"
    # set linestyle
    if cri == "ULS":
        linestyle = "--"  # line style for ULS
    elif cri == "SLS1":
        linestyle = (0, (3, 1, 1, 1))  # line style for SLS1
    elif cri == "SLS2":
        linestyle = ":"  # line style for SLS2
    elif cri == "ENV":
        linestyle = "-"  # line style for ENV
    else:
        linestyle = (0, (1, 10))
    # set linewidth
    if opt == "h":
        linewidth = 0.5
    elif opt == "GWP":
        linewidth = 1.0
    else:
        linewidth = 0.1
    label = sec_typ + ", " + cri + ", optimized for " + opt
#    for j, pl in enumerate(plotted_data):
    for idx, data in enumerate(plotdata):
        plt.subplot(3, 2, idx + 1)
        plt.plot(lengths, data, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
        data_max[idx] = max(data_max[idx], max(data))
for idx, info in enumerate(plotted_data):
    plt.subplot(3, 2, idx + 1)
    plt.xlabel('l [m]')
#    plt.title(info[0])
    plt.ylabel(info[0] + " " + info[1])
    if idx % 2 == 0:
        plt.axis((min(lengths), max(lengths), 0, max(data_max[idx], data_max[idx+1])))
    else:
        plt.axis((min(lengths), max(lengths), 0, max(data_max[idx], data_max[idx-1])))
    plt.legend()
plt.show()


#  VALIDATION
#  isolate cross-sections for verification
validation_idx = 1  # index of length in length-list, corresponding optimal members are separated for further validation
v_members = [member[validation_idx] for member in member_list]

# print some properties of the optimal members, which are useful for manual validation
for idx, member in enumerate(v_members):
    print(legend[idx])
    print("Section Nr. " + str(idx) + " :")
    print(member.section.section_type)
    print("h:")
    print(member.section.h)
    print("admissible load:")
    member.calc_qk_zul_gzt()
    print(member.qk_zul_gzt)
    print("load:")
    print(member.qk)
    print("co2 of section:")
    print(member.section.co2)
    if member.section.section_type == "rc_rec":
        print("x/d:")
        print(member.section.x_p/member.section.d)
        print("di_xu:")
        print(member.section.bw[0])
    print("Admissible deflections (ductile installations):")
    print(member.w_app_adm)
    print("Cracking moment vs. moment")
    print(member.section.mr_p, member.mkd_p)
    print("Calculated deflections (ductile installations):")
    print(member.w_app, member.w_app_ger)
    print(member.a_ed)
    print("fire resistance:")
    member.get_fire_resistance()
    print(member.fire_resistance)
    print(" ")

print("Do manual verification of the data in v_members")
