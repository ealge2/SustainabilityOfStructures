import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import sqlite3  # import modul for SQLite
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Polygon


def plot_dataset(lengths, database_name, criteria, optima, floorstruc, requirements, crsec_type, mat_names,
                 g2k=0.75, qk=2.0, max_iter=50, idx_vrfctn=-1):

    if idx_vrfctn == -1:
        idx_vrfctn = random.randint(0, len(lengths)-1)

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
                section_0 = struct_analysis.RectangularConcrete(concrete, rebar, 1.0, 0.12,
                                                                0.014, 0.15, 0.01, 0.15,
                                                                0.01, 0.15, 2)
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
    vrfctn_members = [[], []]
    for i, members in enumerate(member_list):
        plotdata = [[], [], [], [], []]
        for j, mem in enumerate(members):
            plotdata[0].append(mem.section.h)
            plotdata[1].append(mem.section.h + mem.floorstruc.h)
            plotdata[2].append(mem.section.co2)
            plotdata[3].append(mem.section.co2 + mem.floorstruc.co2)
            plotdata[4].append(mem.section.cost)
            if j == idx_vrfctn:
                vrfctn_members[0].append(mem)
                vrfctn_members[1].append(i)
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
        # plot data
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
            # plot points of verification into graph
            ver_x, ver_y = lengths[idx_vrfctn], data[idx_vrfctn]
            plt.plot(ver_x, ver_y, 'o', color='black', markersize=2)
            plt.annotate(f'#{i}', xy=(ver_x, ver_y),
                         xytext=(ver_x + 0.05*lengths[-1], ver_y),
                         arrowprops=dict(facecolor='black', shrink=0.2, width=0.2, headwidth=2, headlength=4),
                         fontsize=9, color='black', va='center')


    return data_max, vrfctn_members


def plot_section(section):
    # Create a figure and axis
    if section.section_type == "rc_rec":  # Rectangular Reinforced Concrete Cross-Section
        fig, ax = plot_rectangle_with_dimensions(section.b, section.h, 'green', 'x')
        plot_rebars_long(ax, section)
        # add stirrups to plot (if stirrups are defined)
        if section.bw_bg[2]>0:
            plot_stirrups(ax, section)
        legend = (f'{section.concrete_type.mech_prop}, prod_ID:{section.concrete_type.prod_id} \n'
                  f'{section.rebar_type.mech_prop}, prod_ID:{section.rebar_type.prod_id} \n'
                  f'di_xo / s_xo = {section.bw[1][0]:.3f} / {section.bw[1][1]} \n'
                  f'di_xu / s_xu = {section.bw[0][0]:.3f} / {section.bw[0][1]} \n'
                  f'di_bw / s_bw / # = {section.bw_bg[0]} / {section.bw_bg[1]} / {section.bw_bg[2]} \n'
                  f'c_nom = {100*section.c_nom:.1f} cm \n'
                  f'x/d = {section.x_p/section.d:.2f} \n'
                  f'GWP = {section.co2:.0f} kg/m^2')
    elif section.section_type == "wd_rec":  # Rectangular Wooden Cross-Section
        fig, ax = plot_rectangle_with_dimensions(section.b, section.h, 'brown', '/')
        legend = (f'{section.wood_type.mech_prop}, prod_ID:{section.wood_type.prod_id} \n'
                  f'GWP = {section.co2:.0f} kg/m^2')
    else:
        print("no plot for specified section_type defined jet")

    fig.text(0.1, 0.9, legend, ha='left', va='top', fontsize=9, color='black',
             bbox=dict(facecolor='lightgrey', edgecolor='black', boxstyle='round,pad=0.2'))


# Function to plot a rectangular cross-section with given dimensions
def plot_rectangle_with_dimensions(width, height, color='black', hatch='*', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect = patches.Rectangle((offset, offset), width, height, linewidth=1, edgecolor=color, facecolor='none', hatch=hatch,
                             fill=False)

    # Add the rectangle to the plot
    ax.add_patch(rect)

    # Add dimension annotations
    ax.annotate(f'b = {width:.2f} m', xy=(offset + width / 2, 0.05), xytext=(offset + width / 2, 0.06), ha='center')
    ax.annotate(f'h = {height:.2f} m', xy=(0.02, offset + height / 2), xytext=(0.01, offset + height / 2), va='center',
                 rotation='vertical')

    ## Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, width+2*offset)
    ax.set_ylim(0, height+2*offset)

    return fig, ax


def plot_rebars_long(ax, section, color='blue', offset=0.1):
    # get rebar positions
    rebar_positions = get_rebar_positions(section, offset)
    # plot rebars
    for (x, y, r) in rebar_positions:
        rebar = plt.Circle((x, y), r, color='blue')
        ax.add_patch(rebar)


def get_rebar_positions(section, offset):
    # create x and y coordinates of lower longitudinal reinforcement
    y_u = section.h - section.d
    s_xu = section.bw[0][1]
    x_u = [section.b/2]
    while max(x_u) + s_xu < section.b:
        x_u.append(max(x_u) + s_xu)
        x_u.append(min(x_u) - s_xu)

    # create x and y coordinates of upper longitudinal reinforcement
    y_o = section.ds
    s_xo = section.bw[1][1]
    x_o = [section.b/2]
    while max(x_o) + s_xo < section.b:
        x_o.append(max(x_o) + s_xo)
        x_o.append(min(x_o) - s_xo)

    # assemble rebar positions and dimensions relative to cross-section left lower edge corrected by offset
    di_xu = section.bw[0][0]
    di_xo = section.bw[1][0]
    rebar_positions = []
    for xi in x_u:
        rebar_i = (xi+offset, y_u+offset, di_xu)
        rebar_positions.append(rebar_i)
    for xi in x_o:
        rebar_i = (xi+offset, y_o+offset, di_xo)
        rebar_positions.append(rebar_i)

    return rebar_positions

def plot_stirrups(section, offset):
    a = 1