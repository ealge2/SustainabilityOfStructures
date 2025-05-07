# file contains code for generating fundamental plots of first example (simple supported beam)

# IMPORT
import create_dummy_database  # file for creating a "dummy database", as long as no real database is available
import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import matplotlib.pyplot as plt
import plot_datasets  # file with code for plotting results in a standardized way

# INPUT
# create dummy-database
database_name = "dummy_sustainability_1.db"  # define database name
create_dummy_database.create_database(database_name)  # create database

# create materials for reinforced concrete cross-section, derive corresponding design values
# create material for wooden cross-section, derive corresponding design values
timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber1.get_design_values()
timber2 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber2.get_design_values()
timber3 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
timber3.get_design_values()

section_wd0 = struct_analysis.RibWood(timber1, timber2, timber3, 8, 0.12, 0.18, 0.625, 0.027, 0.027)

# create floor structure for solid wooden cross-section
bodenaufbau_holzrippendecke = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                                 ["'Unterlagsboden Zement, 85 mm'", False, False], ["'Glaswolle'", 0.03, False],
                                 ["'Kies gebrochen'", 0.12, False]]
bodenaufbau_wd = struct_analysis.FloorStruc(bodenaufbau_holzrippendecke, database_name)
# create floor structure for solid reinforced concrete cross-section

# define loads on member
g2k = 0.75e3  # n.t. Einbauten
qk = 2e3  # Nutzlast

# define service limit state criteria
req = struct_analysis.Requirements()

# define system lengths for plot
lengths = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

#  define content of plot
to_plot = [[section_wd0, bodenaufbau_wd]]
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
        print(criterion)
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
    elif sec_typ == "rc_rib":
        color = "tab:green"  # color for reinforced concrete
    elif sec_typ == "wd_rib":
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
validation_idx = 6  # index of length in length-list, corresponding optimal members are separated for further validation
v_members = [member[validation_idx] for member in member_list]

# print some properties of the optimal members, which are useful for manual validation
for idx, member in enumerate(v_members):
    print(legend[idx])
    print("Section Nr. " + str(idx) + ":", member.section.section_type)
    print("h:    ", round(member.section.h, 4), 'm')
    print("b:    ", round(member.section.b, 4), 'm')
    print("l:", member.li_max)

    print("admissible load:", member.calc_qk_zul_gzt(), round(member.qk_zul_gzt, 4), 'kN/m2')
    print("load:", round(member.qk,4))
    print("co2 of section:", member.section.co2)
    print("Admissible deflections (ductile installations):")
    print(member.w_app_adm)
    print("Calculated deflections (ductile installations):")
    print(member.w_app)
    print(member.a_ed)
    print("fire resistance:")
    member.get_fire_resistance()
    print(member.fire_resistance)
    print(" ")

print("Do manual verification of the data in v_members")

# # retrieve data from database, find optimal cross-sections and plot results
# mat_names = ["'glue-laminated_timber'", "'solid_structural_timber_(kvh)'"]
#
# data_max, vrfctn_members = plot_datasets.plot_dataset(lengths, database_name, criteria, optima, bodenaufbau_wd, req,
#                                                       "wd_rib", mat_names, g2k, qk, max_iter)
# # plot cross-section of members for verification
# for i, mem in enumerate(vrfctn_members[0]):
#     section = mem.section
#     plot_datasets.plot_section(section)
#     # Show the plot
#     plt.title(f'#{vrfctn_members[1][i]}')
#
# plt.show()
