import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import sqlite3  # import modul for SQLite
import matplotlib.pyplot as plt
from shapely.geometry import Polygon


def plot_dataset(lengths, database_name, criteria, optima, floorstruc, requirements, crsec_type, mat_names,
                 g2k=0.75, qk=2.0, max_iter=50):

    # GENERATE INITIAL CROSS-SECTIONS
    # Search database (table products, attribute material) for products,
    # get prod_id of relevant materials from database and create initial cross-section for each product
    to_plot = []
    connection = sqlite3.connect(database_name)
    cursor = connection.cursor()
    for mat_name in mat_names:
        inquiry = ("SELECT PRO_ID FROM products WHERE"
               " material="+mat_name)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        for i, prod_id in enumerate(result):
            # create materials for wooden cross-sections, derive corresponding design values
            prod_id_str = "'" + str(prod_id[0]) + "'"
            inquiry = ("SELECT mech_prop FROM products WHERE"
                       " PRO_ID=" + prod_id_str)
            cursor.execute(inquiry)
            result = cursor.fetchall()
            mech_prop = "'" + result[0][0] + "'"
            if crsec_type == "wd_rec":
                # create a Wood material object
                timber = struct_analysis.Wood(mech_prop, database_name, prod_id_str)
                timber.get_design_values()
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularWood(timber, 1.0, 0.1, xi=0.02)
            elif crsec_type == "rc_rec":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # create a Rebar material object
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularConcrete(concrete, rebar, 1.0, 0.12, 0.014, 0.15, 0.01, 0.15, 0.01, 0.15)
                # define content of plot-line
            else:
                print("cross-section type is not defined inside function plot_dataset()")
                section_0 = []

            # define content of plot-line
            line_i = [section_0, floorstruc]
            to_plot.append(line_i)

    # ANALYSIS AND OPTIMIZATION OF CROSS-SECTIONS
    member_list = []
    legend = []
    # create plot data
    for i in to_plot:
        for criterion in criteria:
            for optimum in optima:
                members = []
                for length in lengths:
                    sys = struct_analysis.BeamSimpleSup(length)
                    member0 = struct_analysis.Member1D(i[0], sys, i[1], requirements, g2k, qk)
                    opt_section = struct_optimization.get_optimized_section(member0, criterion, optimum, max_iter)
                    opt_member = struct_analysis.Member1D(opt_section, sys, i[1], requirements, g2k, qk)
                    members.append(opt_member)
                member_list.append(members)
                if i[0].section_type[0:2] == "rc":
                    material_lg = i[0].concrete_type.mech_prop + " + " + i[0].rebar_type.mech_prop
                elif i[0].section_type[0:2] == "wd":
                    material_lg = i[0].wood_type.mech_prop
                else:
                    material_lg = "error: section material is not defined"
                legend.append([i[0].section_type, material_lg, criterion, optimum])

    # CREATE DATA OF ENVELOPE AREA OF DATASET
    # create data of envelope area for subplot 1: structural height
    h = [[mem.section.h for mem in sublist] for sublist in member_list]
    h_min = [min(values) for values in zip(*h)]
    h_max = [max(values) for values in zip(*h)]

    # create data of envelope area for subplot 2: total height
    h_tot = [[mem.section.h+mem.floorstruc.h for mem in sublist] for sublist in member_list]
    h_tot_min = [min(values) for values in zip(*h_tot)]
    h_tot_max = [max(values) for values in zip(*h_tot)]

    # create data of envelope area data subplot 3: co2 of structure
    co2 = [[mem.section.co2 for mem in sublist] for sublist in member_list]
    co2_min = [min(values) for values in zip(*co2)]
    co2_max = [max(values) for values in zip(*co2)]

    # create data of envelope area for subplot 4: total co2
    co2_tot = [[mem.section.co2+mem.floorstruc.co2 for mem in sublist] for sublist in member_list]
    co2_tot_min = [min(values) for values in zip(*co2_tot)]
    co2_tot_max = [max(values) for values in zip(*co2_tot)]

    # create data of envelope area for subplot 5: costs of structure
    cost = [[mem.section.cost for mem in sublist] for sublist in member_list]
    cost_min = [min(values) for values in zip(*cost)]
    cost_max = [max(values) for values in zip(*cost)]

    values_min = [h_min, h_tot_min, co2_min, co2_tot_min, cost_min]
    values_max = [h_max, h_tot_max, co2_max, co2_tot_max, cost_max]

    # PLOT DATASET TO FIGURE
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
        sec_typ, mat, cri, opt = legend[i]
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
        label = sec_typ + ", " + mat + ", " + cri + ", optimized for " + opt
        for idx, data in enumerate(plotdata):
            plt.subplot(3, 2, idx + 1)
            # prepare area
            coords = list(zip(lengths, values_max[idx])) + list(zip(lengths[::-1], values_min[idx][::-1]))
            # create a polygon from the coordinates
            polygon = Polygon(coords)
            # extract the x and y coordinates for plotting
            x, y = polygon.exterior.xy
            # plot area
            plt.fill(x, y, alpha=0.2, facecolor=color)
            # plot lines
            plt.plot(lengths, data, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
            data_max[idx] = max(data_max[idx], max(data))

    return data_max