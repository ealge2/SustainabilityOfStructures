import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import sqlite3  # import modul for SQLite
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Polygon
from scipy.interpolate import interp1d
from scipy.spatial import ConvexHull



# PLOT DATASETS OF MEMBERS WITH DEFINED CROSS_SECTIONS AND VARIED MATERIALS
# ----------------------------------------------------------------------------------------------------------------------
def plot_dataset(lengths, database_name, criteria, optima, floorstruc, requirements, crsec_type, mat_names,
                 g2k=0.75, qk=2.0, max_iter=100, idx_vrfctn=-1):

    if idx_vrfctn == -1:
        idx_vrfctn = random.randint(0, len(lengths)-1)

    # GENERATE INITIAL CROSS-SECTIONS
    # Search database (table products, attribute material) for products,
    # get prod_id of relevant materials from database and create initial cross-section for each product
    to_plot = []
    connection = sqlite3.connect(database_name)
    cursor = connection.cursor()
    for mat_name in mat_names:
        inquiry = ("""
                SELECT PRO_ID FROM products
                WHERE density IS NOT NULL
                AND mech_prop IS NOT NULL
                AND "material [string]" LIKE """ + mat_name
        )
        # inquiry = ("SELECT PRO_ID FROM products WHERE"
        #            " material=" + mat_name)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        for i, prod_id in enumerate(result):
            # create materials for wooden cross-sections, derive corresponding design values
            prod_id_str = "'" + str(prod_id[0]) + "'"
            inquiry = ("""
                    SELECT mech_prop FROM products
                    WHERE  PRO_ID LIKE """ + prod_id_str
            )
            # inquiry = ("SELECT mech_prop FROM products WHERE"
            #            " PRO_ID=" + prod_id_str)
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
                # ToDo integrate calculations with low and high value for rebar emissions
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularConcrete(concrete, rebar, 1.0, 0.12,
                                                                0.014, 0.15, 0.01, 0.15,
                                                                0, 0.15, 2)
                color = "tab:green"


            elif crsec_type == "rc_rib":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # create a Rebar material object
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                #ToDo integrate calculations with low and high value for rebar emissions

                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RibbedConcrete(concrete, rebar, 4, 1.0, 0.14, 0.3, 0.18, 0.01, 0.15, 0.01, 0.15, 0.02, 2, 0.01, 0.15, 2)

            elif crsec_type == "wd_rib":
                # create a Wood material object
                timber1 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
                timber1.get_design_values()
                timber2 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
                timber2.get_design_values()
                timber3 = struct_analysis.Wood("'GL24h'", database_name)  # create a Wood material object
                timber3.get_design_values()
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RibWood(timber1, timber2, timber3, 4, 0.12, 0.18, 0.625, 0.027, 0.027)



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
                elif i[0].section_type == "wd":
                    material_lg = i[0].wood_type.mech_prop
                elif i[0].section_type == "wd_rib":
                    material_lg = i[0].wood_type_1.mech_prop
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
        elif sec_typ == "rc_rib":
            color = "tab:green"  # color for wood
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

# PLOT GEOMETRY OF SECTIONS
# ----------------------------------------------------------------------------------------------------------------------
def plot_section(section):
    # Create a figure and axis
    if section.section_type == "rc_rec":  # Rectangular Reinforced Concrete Cross-Section
        fig, ax, offset = plot_rectangle_with_dimensions(section.b, section.h, 'green', 'x')
        plot_rebars_long(ax, section, offset)
        # add stirrups to plot (if stirrups are defined)
        if section.bw_bg[0] > 0 and section.bw_bg[2] > 0:
            plot_stirrups(ax, section, offset)
        legend = (f'{section.concrete_type.mech_prop}, prod_ID:{section.concrete_type.prod_id} \n'
                  f'{section.rebar_type.mech_prop}, prod_ID:{section.rebar_type.prod_id} \n'
                  f'di_xo / s_xo = {section.bw[1][0]:.3f} / {section.bw[1][1]} \n'
                  f'di_xu / s_xu = {section.bw[0][0]:.3f} / {section.bw[0][1]} \n'
                  f'di_stir / s_stir / n = {section.bw_bg[0]} / {section.bw_bg[1]} / {section.bw_bg[2] }\n'
                  f'c_nom = {100*section.c_nom:.1f} cm \n'
                  f'x/d = {section.x_p/section.d:.2f} \n'
                  f'GWP = {section.co2:.0f} kg/m^2')
    elif section.section_type == "wd_rec":  # Rectangular Wooden Cross-Section
        fig, ax, offset = plot_rectangle_with_dimensions(section.b, section.h, 'brown', '/')
        legend = (f'{section.wood_type.mech_prop}, prod_ID:{section.wood_type.prod_id} \n'
                  f'GWP = {section.co2:.0f} kg/m^2')

    elif section.section_type == "rc_rib": #Betonrippenquerschnitte
        fig, ax, offset = plot_rib_with_dimensions(section.b, section.b_w, section.h, section.h_f, 'green', 'x')
        legend = (f'{section.concrete_type.mech_prop}, prod_ID:{section.concrete_type.prod_id} \n'
              f'length = {section.l0} \n'
              f'{section.rebar_type.mech_prop}, prod_ID:{section.rebar_type.prod_id} \n'
              f'di_r = {section.bw_r[0]:.3f} \n'
              f'di_xu / s_xu = {section.bw[0][0]:.3f} / {section.bw[0][1]} \n'
              #f'di_stir / s_stir / n = {section.bw_bg[0]} / {section.bw_bg[1]} / {section.bw_bg[2]}\n'
              f'c_nom = {100 * section.c_nom:.1f} cm \n'
              f'x/d = {section.x_p / section.d:.2f} \n'
              f'h, hf, hw, b, bw = {section.h:.2f}, {section.h_f:.2f}, {section.h_w:.2f}, {section.b:.2f}, {section.b_w:.2f} \n'
              f'GWP = {section.co2:.0f} kg/m^2')

    elif section.section_type == "wd_rib": #Betonrippenquerschnitte
        fig, ax, offset = plot_wd_rib_with_dimensions(section.b, section.h, section.a, section.t2, section.t3, 'brown', 'x')
        legend = (f'{section.wood_type_1.mech_prop}, prod_ID:{section.wood_type_1.prod_id} \n'
              f'length = {section.l0} \n'
              f'h, b, a, t2, t3 = {section.h:.2f}, {section.b:.2f}, {section.a:.2f}, {section.t2:.2f}, {section.t3:.2f} \n'
              f'GWP = {section.co2:.0f} kg/m^2')

    else:
        print("no plot for specified section_type defined jet")
        fig, ax = plt.subplots()
        legend = f'no plot for section_type "{section.section_type}" defined jet'
    fig.text(0.01, 0.99, legend, ha='left', va='top', fontsize=9, color='black',
             bbox=dict(facecolor='lightgrey', edgecolor='black', boxstyle='round,pad=0.2'))


# Function to plot a rectangular cross-section with given dimensions
def plot_rectangle_with_dimensions(width, height, color='black', hatch='*', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect = patches.Rectangle((offset, offset), width, height, linewidth=1, edgecolor=color, facecolor='none',
                             hatch=hatch, fill=False)

    # Add the rectangle to the plot
    ax.add_patch(rect)

    # Add dimension annotations
    ax.annotate(f'b = {width:.2f} m', xy=(offset + width / 2, 0.05), xytext=(offset + width / 2, 0.06), ha='center')
    ax.annotate(f'h = {height:.2f} m', xy=(0.02, offset + height / 2), xytext=(0.01, offset + height / 2),
                va='center', rotation='vertical')

    # Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, width+5*offset)
    ax.set_ylim(0, height+4*offset)

    return fig, ax, offset

def plot_rib_with_dimensions(b, bw, h, hf, color='black', hatch='*', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect_flange = patches.Rectangle((offset, offset+h-hf), 2*b, hf, linewidth=1, edgecolor=color, facecolor='none',
                             hatch=hatch, fill=False)
    rect_rib1 = patches.Rectangle((offset+b/2-bw/2, offset), bw, h-hf, linewidth=1, edgecolor=color, facecolor='none',
                             hatch=hatch, fill=False)
    rect_rib2 = patches.Rectangle((offset + 3*b / 2 - bw / 2, offset), bw, h - hf, linewidth=1, edgecolor=color,
                                 facecolor='none',
                                 hatch=hatch, fill=False)

    # Add the rectangle to the plot
    ax.add_patch(rect_flange)
    ax.add_patch(rect_rib1)
    ax.add_patch(rect_rib2)

    # # Add dimension annotations
    # ax.annotate(f'b = {width:.2f} m', xy=(offset + width / 2, 0.05), xytext=(offset + width / 2, 0.06), ha='center')
    # ax.annotate(f'h = {height:.2f} m', xy=(0.02, offset + height / 2), xytext=(0.01, offset + height / 2),
    #             va='center', rotation='vertical')

    # Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, b + 10 * offset)
    ax.set_ylim(0, h + 4 * offset)

    return fig, ax, offset

def plot_wd_rib_with_dimensions(b, h, a, t2, t3, color='black', hatch='--', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect_flange2 = patches.Rectangle((offset, offset), 2*a, t2, linewidth=1, edgecolor=color, facecolor='none',
                             hatch='--', fill=False)
    rect_flange3 = patches.Rectangle((offset, offset + t2+h), 2*a, t3, linewidth=1, edgecolor=color, facecolor='none',
                                     hatch='--', fill=False)
    rect_rib1 = patches.Rectangle((offset+a/2, offset+t2), b, h, linewidth=1, edgecolor=color, facecolor='none',
                             hatch='-', fill=False)
    rect_rib2 = patches.Rectangle((offset+3*a/2, offset+t2), b, h, linewidth=1, edgecolor=color, facecolor='none',
                             hatch='-', fill=False)


    # Add the rectangle to the plot
    ax.add_patch(rect_flange2)
    ax.add_patch(rect_flange3)
    ax.add_patch(rect_rib1)
    ax.add_patch(rect_rib2)

    # Add dimension annotations
    ax.annotate(f'b_Rippe = {b:.2f} m', xy=(offset + b, 0.05), xytext=(offset + b / 2, 0.06), ha='center')
    ax.annotate(f'b = {b:.2f} m', xy=(offset + b/2, 0.05), xytext=(offset + 3*b / 2, -0.16), ha='center')
    ax.annotate(f'h = {h:.2f} m', xy=(0.02, offset + h / 2), xytext=(0.01, offset + h / 2),
                va='center', rotation='vertical')

    # #Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, 2*a + 2 * offset)
    ax.set_ylim(0, h+t2+t3 + 4 * offset)

    return fig, ax, offset

def plot_rebars_long(ax, section, offset, color='blue'):
    # get rebar positions
    rebar_positions = get_rebar_positions(section)
    # plot rebars
    for (x, y, di) in rebar_positions:
        rebar = plt.Circle((x+offset, y+offset), di/2, color=color)
        ax.add_patch(rebar)

def plot_stirrups(ax, section, offset, color='blue'):
    # get stirrup positions
    stirrup_positions = get_stirrup_positions(section)
    # plot stirrups
    for (x1, y1, x2, y2) in stirrup_positions:
        stirrup = plt.Line2D([x1+offset, x2+offset], [y1+offset, y2+offset], color=color, linewidth=2)
        ax.add_line(stirrup)

def get_rebar_positions(section):
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

    # assemble rebar positions and dimensions relative to cross-section left lower edge
    di_xu = section.bw[0][0]
    di_xo = section.bw[1][0]
    rebar_positions = []
    for xi in x_u:
        rebar_i = (xi, y_u, di_xu)
        rebar_positions.append(rebar_i)
    for xi in x_o:
        rebar_i = (xi, y_o, di_xo)
        rebar_positions.append(rebar_i)
    return rebar_positions


def get_stirrup_positions(section):
    # create coordinates of stirrup edge points
    n_stirrup = section.bw_bg[2]
    di_stirrup = section.bw_bg[0]
    edge_dist = section.c_nom + di_stirrup/2
    # create x-coordinates of vertical stirrups
    x_linspace = np.linspace(edge_dist, section.b-edge_dist, n_stirrup)
    # create y-coordinates of stirrups-endings
    y_u = edge_dist
    y_o = section.h - edge_dist
    # create lines between edge points
    stirrup_positions = []
    for idx, xi in enumerate(x_linspace):
        # create vertical lines
        line_vert = (xi, y_u, xi, y_o)
        stirrup_positions.append(line_vert)
        if idx < x_linspace.size-1:
            # create horizontal lines
            line_u = (xi, y_u, x_linspace[idx+1], y_u)
            stirrup_positions.append(line_u)
            line_o = (xi, y_o, x_linspace[idx+1], y_o)
            stirrup_positions.append(line_o)
    return stirrup_positions


# PLOT DATASETS OF CROSS_SECTION WITH VARIED MATERIALS IN M-CHI PLOT AND PLOT THE OPTIMIZED SECTIONS FOR VALIDATION
# ----------------------------------------------------------------------------------------------------------------------
def plot_section_dataset(database_name, crsec_type, mat_names, ax, gwp_budget=50):
    # GENERATE INITIAL CROSS-SECTIONS
    # Search database (table products, attribute material) for products,
    # get prod_id of relevant materials from database and create initial cross-section for each product
    to_plot = []
    connection = sqlite3.connect(database_name)
    cursor = connection.cursor()
    mat_nr = -1
    x_values = []
    y_values = []
    for mat_name in mat_names:
        inquiry = ("SELECT PRO_ID FROM products WHERE"
                   " material=" + mat_name)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        for i, prod_id in enumerate(result):
            mat_nr += 1  # number for annotation in plot
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
                section_0 = struct_analysis.RectangularWood(timber, 1.0, 0.12)
                color = "tab:brown"
            elif crsec_type == "rc_rec":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # create a Rebar material object
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularConcrete(concrete, rebar, 1.0, 0.12,
                                                                0.014, 0.15, 0.01, 0.15,
                                                                0, 0.15, 2)
                color = "tab:green"
            elif crsec_type == "rc_rib":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # create a Rebar material object
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RibbedConcrete(concrete, rebar, 4, 1.0, 0.14, 0.3, 0.18, 0.01, 0.15, 0.01, 0.15, 0.02, 2, 0.01, 0.15, 2)
                color = "tab:green"
## XXXXXXXXXXX neuen Querschnittstyp initialisieren
            else:
                print("cross-section type is not defined inside function plot_dataset()")
                section_0 = []
                color = "tab:grey"


            #
            # maximizing Mu by varying the geometry within the max allowed gwp-budget.
            opt_section = struct_optimization.get_opt_sec(section_0, gwp_budget)

            # add M-Chi relationship to plot
            x_values, y_values = plot_m_chi(opt_section, ax, mat_nr, x_values, y_values)

            # plot cross-section
            plot_section(opt_section)
            plt.title(f'#{mat_nr}')

    # ## Add envelope of dataset to plot
    # # Combine x and y into a single array for ConvexHull
    # x = [item for sublist in x_values for item in sublist]
    # y = [item for sublist in y_values for item in sublist]
    # points = np.column_stack((x, y))
    # # points = [item for sublist in points_nested for item in sublist]
    # print(points.shape)
    # hull = ConvexHull(points)
    #
    # # Plot the envelope area (convex hull)
    # for simplex in hull.simplices:
    #     ax.plot(points[simplex, 0], points[simplex, 1], 'r-')
#
#     # Fill the hull area
#     hull_vertices = points[hull.vertices]
#     ax.fill(hull_vertices[:, 0], hull_vertices[:, 1], color, alpha=0.3, label="Envelope")

    ## Add envelope of dataset to plot
    # Define a common x-axis for interpolation
    common_x = sorted(set(np.concatenate(x_values)))
    # Interpolate y-values onto the common x-axis
    interpolated_y_values_min = []
    interpolated_y_values_max = []
    # y_test = []
    for x, y in zip(x_values, y_values):
        # define interpolation functions with different rules, when x is out of range
        f1 = interp1d(x, y, kind='linear', bounds_error=False, fill_value="extrapolate")
        f2 = interp1d(x, y, kind='linear', bounds_error=False, fill_value=(0, 0))
        a = []
        b = []
        for xi in common_x:
            if xi >= 0:  # apply interpolation functions withe rules for positive x_values
                a.append(float(f1(xi)))
                b.append(float(f2(xi)))
            else:  # apply interpolation functions withe rules for negative x_values
                b.append(float(f1(xi)))
                a.append(float(f2(xi)))
        interpolated_y_values_min.append(a)
        interpolated_y_values_max.append(b)

    y_lower = np.min(interpolated_y_values_min, axis=0)
    y_upper = np.max(interpolated_y_values_max, axis=0)

    # Fill the area between the envelope
    ax.fill_between(common_x, y_lower, y_upper, color=color, alpha=0.3, label='Envelope')



# plot m-chi relationship for a defined cross-section
def plot_m_chi(section, ax, i, x_values, y_values):
    # add M-Chi relationship to plot
    if section.section_type[0:2] == "wd":
        color = "tab:brown"  # color for wood
        chi_u_p, chi_u_n = section.mu_max/section.ei1, section.mu_min/section.ei1
        x = [chi_u_n, 0, chi_u_p]
        y = [section.mu_min, 0, section.mu_max]
    elif section.section_type[0:2] == "rc":
        color = "tab:green"  # color for reinforced concrete
        chi_r1_p, chi_r1_n = section.mr_p / section.ei1, section.mr_n / section.ei1
        chi_r2_p, chi_r2_n = section.mr_p / section.ei2, section.mr_n / section.ei2
        chi_y_p, chi_y_n = section.mu_max / section.ei2, section.mu_min / section.ei2
        chi_u_p, chi_u_n = section.concrete_type.ec2d/section.x_p, -section.concrete_type.ec2d/section.x_n
        if section.mr_p <= section.mu_max and section.mr_n >= section.mu_min:  # ductile behavior in both directions
            x = [chi_u_n, chi_y_n, chi_r2_n, chi_r1_n, 0, chi_r1_p, chi_r2_p, chi_y_p, chi_u_p]
            y = [section.mu_min, section.mu_min, section.mr_n, section.mr_n, 0, section.mr_p, section.mr_p,
                 section.mu_max, section.mu_max]
        elif section.mr_p <= section.mu_max and section.mr_n < section.mu_min:  # ductile behavior only for pos. moments
            x = [chi_r1_n, chi_r1_n, 0, chi_r1_p, chi_r2_p, chi_y_p, chi_u_p]
            y = [0, section.mr_n, 0, section.mr_p, section.mr_p,
                 section.mu_max, section.mu_max]
        elif section.mr_p > section.mu_max and section.mr_n >= section.mu_min:  # ductile behavior only for neg. moments
            x = [chi_u_n, chi_y_n, chi_r2_n, chi_r1_n, 0, chi_r1_p, chi_r1_p]
            y = [section.mu_min, section.mu_min, section.mr_n, section.mr_n, 0, section.mr_p, 0]
        else:  # no ductile behavior
            x = [chi_r1_n, chi_r1_n, 0, chi_r1_p, chi_r1_p]
            y = [0, section.mr_n, 0, section.mr_p, 0]
    else:
        print("M-Chi plot is not defined for sections of type " + section.section_type)
        x, y = 0, 0
        color = "tab:black"
    y_mod = [yi/1e3 for yi in y]  # modify unit of moment from Nm/m to kNm/m
    ax.plot(x, y_mod, color=color)  # plot m-Chi relationship with unit of moment: kNm/m
    ax.annotate(f'#{i}', xy=(x[-1], y_mod[-1]), xytext=(x[-1]*1.1, y_mod[-1]),
                 arrowprops=dict(facecolor='black', shrink=0.2, width=0.2, headwidth=2, headlength=4),
                 fontsize=9, color='black', va='center')
    x_values.append(x)
    y_values.append(y_mod)
    return x_values, y_values
    ## Define and plot envelope


# plt.annotate(f'#{i}', xy=(ver_x, ver_y),
#                          xytext=(ver_x + 0.05*lengths[-1], ver_y),
#                          arrowprops=dict(facecolor='black', shrink=0.2, width=0.2, headwidth=2, headlength=4),
#                          fontsize=9, color='black', va='center')