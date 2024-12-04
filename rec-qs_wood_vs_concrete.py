# file contains code for generating plots for "Austausch Nr. 1 mit erweitertem Projektteam"
# parts of the code will to be integrated in a file, which contains general code for automated plot generation

import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import matplotlib.pyplot as plt

# max. number of iterations per optimization. Fast results: max_iterations = 50, good results: max iterations = 1000
# max_iter = 1000

# create dummy-database
database_name = "dummy_sustainability.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# create material for wooden cross-section, derive corresponding design values
timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber1.get_design_values()
# create materials for reinforced concrete cross-section, derive corresponding design values
concrete1 = struct_analysis.ReadyMixedConcrete("'C25/30'", database_name)
concrete1.get_design_values()
reinfsteel1 = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
reinfsteel1.get_design_values()

# create initial wooden rectangular cross-section
section_wd0 = struct_analysis.RectangularWood(timber1, 1.0, 0.1)
# create initial reinforced concrete rectangular cross-section
section_rc0 = struct_analysis.RectangularConcrete(concrete1, reinfsteel1, 1.0, 0.1, 0.012, 0.15, 0.01, 0.15)

# create floor structure for solid wooden cross-section
bodenaufbau_brettstappeldecke = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                                 ["'Unterlagsboden Zement, 85 mm'", False, False], ["'Glaswolle'", 0.03, False],
                                 ["'Kies gebrochen'", 0.12, False]]
bodenaufbau_wd = struct_analysis.FloorStruc(bodenaufbau_brettstappeldecke, database_name)
# create floor structure for solid reinforced concrete cross-section
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
lengths = [3, 4, 5, 6, 7, 8, 9, 10, 12]

#  define content of plot
to_plot = [[section_rc0, bodenaufbau_rc], [section_wd0, bodenaufbau_wd]]
criteria = ["ULS", "SLS1"]
optima = ["GWP", "h"]
plotted_data = [["h_struct", "[m]"], ["h_tot", "[m]"], ["GWP_struct", "[kg-CO2-eq]"], ["GWP_tot", "[kg-CO2-eq]"],
                ["cost_struct", "[CHF]"]]
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
                opt_section = struct_optimization.get_optimized_section(member0, criterion, optimum)
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
        plotdata[4].append(mem.section.cost + mem.floorstruc.co2)
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
        linestyle = "-"  # line style for ULS
    elif cri == "SLS1":
        linestyle = "--"  # line style for SLS1
    else:
        linestyle = ""
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

#isolate cross-sections for verification
idx_length = 5
v_members = [member[5] for member in member_list]

for idx, member in enumerate(v_members):
    print(legend[idx])
    print("Section Nr. " + str(idx) + " :")
    print(member.section.section_type)
    print("admissible load:")
    member.calc_qk_zul_gzt()
    print(member.qk_zul_gzt)
    print("x/d:")
    print(member.section.x_p/member.section.d)
    print("Admissible deflections (ductile installations):")
    print(member.w_install_adm)
    print("Calculated deflections (ductile installations):")
    print(member.w_install)


#  XXXXXXXXXXXXXXXXXXXXXXX Optimization has to be fixed for some cross-sectoins!XXXXXXXXXXXXXXX
print("Do manual verification of the data in v_members")

print("End of code")
