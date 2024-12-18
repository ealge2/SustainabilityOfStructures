# File enthält Code für die Strukturanalyse (Bauteil- und Querschnittsanalyse)
# units: [m], [kg], [s], [N], [CHF]

# Abgebildete Materialien:
# - Beton
# - Betonstahl
# - Holz
#
# Abgebildete Querschnitte 1D:
# - Betonrechteck-QS
# - Holzrechteck-QS
#
# Abgebildete Statische Systeme 1D:
# - Einfacher Balken
#
# Weitere Klassen:
# - Bauteil 1D
# - Bodenaufbauschicht
# - Bodenaufbau
# - Rechteckquerschnitte

import sqlite3  # import modul for SQLite
import numpy as np

class Wood:
    # defines properties of wooden material
    def __init__(self, mech_prop, database):  # retrieve basic mechanical data from database
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = ("SELECT strength_bend, strength_shea, E_modulus, density_load FROM material_prop WHERE"
                   " name="+mech_prop)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fmk, self.fvd, self.Emmean, self.weight = result[0]
        # get GWP properties from database
        inquiry = "SELECT density, GWP, cost, cost2 FROM products WHERE mech_prop="+mech_prop
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.density, self.GWP, self.cost, self.cost2 = result[0]
        self.fmd = float()

    def get_design_values(self, gamma_m=1.7, eta_m=1, eta_t=1, eta_w=1):  # calculate design values
        if self.mech_prop[1:3] == "GL":
            gamma_m = 1.5  # SIA 265, 2.2.5: reduzierter Sicherheitsbeiwert für BSH

        self.fmd = self.fmk * eta_m * eta_t * eta_w / gamma_m  # SIA 265, 2.2.2, Formel (3)


class ReadyMixedConcrete:
    # defines properties of concrete material
    def __init__(self, mech_prop, database):  # retrieve basic mechanical data from database (self, table,
        self.ec2d = float()
        self.tcd = float()
        self.fcd = float()
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = ("SELECT strength_comp, strength_tens, E_modulus, density_load FROM material_prop WHERE name="
                   + mech_prop)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fck, self.fctm, self.Ecm, self.weight = result[0]
        # get GWP properties from database
        inquiry = "SELECT density, GWP, cost, cost2 FROM products WHERE mech_prop="+mech_prop
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.density, self.GWP, self.cost, self.cost2 = result[0]

    def get_design_values(self, gamma_c=1.5, eta_t=1):  # calculate design values
        eta_fc = min((30e6/self.fck) ** (1/3), 1)  # SIA 262, 4.2.1.2, Formel (26)
        self.fcd = self.fck * eta_fc * eta_t / gamma_c  # SIA 262, 2.3.2.3, Formel (2)
        self.tcd = 0.3 * eta_t * self.fck ** 0.5/gamma_c  # SIA 262, 2.3.2.4, Formel (3)
        self.ec2d = 0.003  # SIA 262, 4.2.4, Tabelle 8


class SteelReinforcingBar:
    # defines properties of reinforcement  material
    def __init__(self, mech_prop, database):
        # retrieve basic mechanical data from database (self, table, database name)
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = "SELECT strength_tens, E_modulus FROM material_prop WHERE name="+mech_prop
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fsk, self.Es = result[0]
        # get GWP properties from database
        inquiry = "SELECT density, GWP, cost FROM products WHERE mech_prop="+mech_prop
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.density, self.GWP, self.cost = result[0]
        self.fsd = float()

    def get_design_values(self, gamma_s=1.15):  # calculate design values
        self.fsd = self.fsk/gamma_s  # SIA 262, 2.3.2.5, Formel (4)


class Section:
    # contains fundamental section properties like sectoion type weight, resistance and stiffness
    def __init__(self, section_type):
        self.section_type = section_type
        # self.mu_max = float
        # self.mu_min = float
        # self.vu = float
        # self.qs_class_n = int
        # self.qs_class_p = int
        # self.g0k = float
        # self.ei1 = float
        # self.co2 = float
        # self.cost = float


class SupStrucRectangular(Section):
    # defines cross-section dimensions and has methods to calculate static properties of rectangular,
    # non-cracked sections

    def __init__(self, section_type, b, h, phi=0):  # create a rectangular object
        super().__init__(section_type)
        self.b = b  # width [m]
        self.h = h  # height [m]
        self.a_brutt = self.calc_area()
        self.iy = self.calc_moment_of_inertia()
        self.phi = phi

    def calc_area(self):
        #  in: width b [m], height h [m]
        #  out: area [m^2]
        a_brutt = self.b * self.h
        return a_brutt

    def calc_moment_of_inertia(self):
        #  in: width b [m], height h [m]
        #  out: second moment of inertia Iy [m^4]
        iy = self.b * self.h ** 3 / 12
        return iy

    def calc_strength_elast(self, fy, ty):
        #  in: yielding strength fy [Pa], shear strength ty [Pa]
        #  out: elastic bending resistance [Nm], elastic shear resistance [N]
        mu_el = self.iy * fy * 2 / self.h
        vu_el = self.b * self.h * ty / 1.5
        return mu_el, vu_el

    def calc_strength_plast(self, fy, ty):
        #  in: yielding strength fy [Pa], shear strength ty [Pa]
        #  out: plastic bending resistance [Nm], plastic shear resistance [N]
        mu_pl = self.b * self.h ** 2 * fy / 4
        vu_pl = self.b * self.h * ty
        return mu_pl, vu_pl

    def calc_weight(self, spec_weight):
        #  in: specific weight [N/m^3]
        #  out: weight of cross section per m length [N/m]
        w = spec_weight * self.a_brutt
        return w

#Ripped cross sections
class SupStrucRipped(Section):




class RectangularWood(SupStrucRectangular, Section):
    # defines properties of rectangular, wooden cross-section
    def __init__(self, wood_type, b, h, phi=0.6):  # create a rectangular timber object
        section_type = "wd_rec"
        super().__init__(section_type, b, h, phi)
        self.wood_type = wood_type
        mu_el, vu_el = self.calc_strength_elast(wood_type.fmd, wood_type.fvd)
        self.mu_max, self.mu_min = [mu_el, mu_el]   #Readme: Why is this needed for wood?
        self.vu = vu_el
        self.qs_class_n, self.qs_class_p = [3, 3]   #Readme: What is this needed for?
        self.g0k = self.calc_weight(wood_type.weight)
        self.ei1 = self.wood_type.Emmean*self.iy  # elastic stiffness [Nm^2]
        self.co2 = self.a_brutt * self.wood_type.GWP * self.wood_type.density  # [kg_CO2_eq/m]
        self.cost = self.a_brutt * self.wood_type.cost


class RectangularConcrete(SupStrucRectangular):
    # defines properties of rectangular, reinforced concrete cross-section
    def __init__(self, concrete_type, rebar_type, b, h, di_xu, s_xu, di_xo, s_xo, phi=2.0, c_nom=0.03):

        # create a rectangular concrete object
        section_type = "rc_rec"
        super().__init__(section_type, b, h, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom
        self.bw = [[di_xu, s_xu], [di_xo, s_xo]]
        # self.bw_bg = XXXXXXXXXXToDoXXXXXXXXXX
        [self.d, self.ds] = self.calc_d()
        [self.mu_max, self.x_p, self.as_p, self.qs_class_p] = self.calc_mu('pos')
        [self.mu_min, self.x_n, self.as_n, self.qs_class_n] = self.calc_mu('neg')
        # [self.vu, self.as_bg] = self.calc_shear_resistance() XXXXXXXXXXToDoXXXXXXXXXX
        self.g0k = self.calc_weight(concrete_type.weight)
        a_s_tot = self.as_p + self.as_n  # add area of stirrups XXXXXXXXXXToDoXXXXXXXXXX
        co2_rebar = a_s_tot * self.rebar_type.GWP * self.rebar_type.density  # [kg_CO2_eq/m]
        co2_concrete = (self.a_brutt-a_s_tot) * self.concrete_type.GWP * self.concrete_type.density  # [kg_CO2_eq/m]
        self.ei1 = self.concrete_type.Ecm*self.iy  # elastic stiffness concrete (uncracked behaviour) [Nm^2]
        self.co2 = co2_rebar + co2_concrete
        self.cost = (a_s_tot * self.rebar_type.cost + (self.a_brutt-a_s_tot) * self.concrete_type.cost
                     + self.concrete_type.cost2)
    #   self.ei2 = # XXXXXXXXXXToDoXXXXXXXXXX

    def calc_d(self):
        d = self.h - self.c_nom - self.bw[0][0]/2
        ds = self.h - self.c_nom - self.bw[1][0]/2
        return d, ds

    def calc_mu(self, sign='pos'):
        b = self.b
        fsd = self.rebar_type.fsd
        fcd = self.concrete_type.fcd
        if sign == 'pos':
            [mu, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[0][0], self.bw[0][1], self.d, b, fsd, fcd)
        elif sign == 'neg':
            [mu, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[1][0], self.bw[1][1], self.ds, b, fsd, fcd)
        else:
            [mu, x, a_s, qs_klasse] = [0, 0, 0, 0]
            print("sigen of moment resistance has to be 'neg' or 'pos'")
        return mu, x, a_s, qs_klasse

    @staticmethod
    def mu_unsigned(di, s, d, b, fsd, fcd):
        # units input: [m, m, m, m, N/m^2, N/m^2]
        a_s = np.pi * di ** 2 / (4 * s) * b  # [m^2]
        omega = a_s * fsd / (d * b * fcd)  # [-]
        mu = a_s * fsd * d * (1-omega/2)  # [Nm]
        x = omega * d / 0.85  # [m]
        if x/d <= 0.35:
            return mu, x, a_s, 1
        elif x/d <= 0.5:
            return mu, x, a_s, 2
        else:
            return mu, x, a_s, 99  # Querschnitt hat ungenügendes Verformungsvermögen


class MatLayer:  # create a material layer
    def __init__(self, mat_name, h_input, roh_input, database):  # get initial data from database
        self.name = mat_name
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get properties from database
        inquiry = "SELECT h_fix, density, weight, GWP FROM floor_struc_prop WHERE name=" + mat_name
        cursor.execute(inquiry)
        result = cursor.fetchall()
        h_fix, density, weight, self.GWP = result[0]
        if h_input is False:
            self.h = h_fix
        else:
            self.h = h_input
        if roh_input is False:
            self.density = density
            self.weight = weight
        else:
            self.density = roh_input
            self.weight = roh_input*10
        self.gk = self.weight * self.h  # weight per area in N/m^2
        self.co2 = self.density * self.h * self.GWP  # CO2-eq per area in kg-C02/m^2


class FloorStruc:  # create a floor structure
    def __init__(self, mat_layers, database_name):
        self.layers = []
        self.co2 = 0
        self.gk_area = 0
        self.h = 0
        for mat_name, h_input, roh_input in mat_layers:
            current_layer = MatLayer(mat_name, h_input, roh_input, database_name)
            self.layers.append(current_layer)
            self.co2 += current_layer.co2
            self.gk_area += current_layer.gk
            self.h += current_layer.h


class BeamSimpleSup:
    def __init__(self, length):
        self.l_tot = length
        self.li_max = self.l_tot  # max span (used for calculation of admissible deflections)
        self.alpha_m = [0, 1/8]
        self.qs_cl_erf = [3, 3]  # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 5/384


class Member1D:
    def __init__(self, section, system, floorstruc, requirements, g2k=0.0, qk=2e3, psi0=0.7, psi1=0.5, psi2=0.3):
        self.section = section
        self.system = system
        self.floorstruc = floorstruc
        self.requirements = requirements
        self.g0k = self.section.g0k
        self.g1k = self.floorstruc.gk_area
        self.g2k = g2k
        self.gk = self.g0k + self.g1k + self.g2k
        self.qk = qk
        self.psi = [psi0, psi1, psi2]
        self.q_rare = self.gk + self.qk
        self.q_freq = self.gk + self.psi[1]*self.qk
        self.q_per = self.gk + self.psi[2]*self.qk
        self.w_install_adm = self.system.li_max/self.requirements.lw_install
        self.w_use_adm = self.system.li_max/self.requirements.lw_use
        self.w_app_adm = self.system.li_max/self.requirements.lw_app
        self.qu = self.calc_qu()
        self.qk_zul_gzt = float

        # calculation of deflections (uncracked cross-section, method for cracked cross-section is not implemented jet)
        if self.requirements.install == "ductile":
            self.w_install = self.system.alpha_w * (
                        self.q_freq + self.q_per * (self.section.phi - 1)) * self.system.l_tot ** 4 / self.section.ei1
        elif self.requirements.install == "brittle":
            self.w_install = self.system.alpha_w * (
                    self.q_rare + self.q_per * (self.section.phi - 1)) * self.system.l_tot ** 4 / self.section.ei1
        self.w_use = self.system.alpha_w * (
                    self.q_freq - self.gk) * self.system.l_tot ** 4 / self.section.ei1
        self.w_app = self.system.alpha_w * (
                self.q_per * (1 + self.section.phi)) * self.system.l_tot ** 4 / self.section.ei1
        self.co2 = system.l_tot * (floorstruc.co2 + section.co2)

    def calc_qu(self):
        # calculates maximal load qu in respect to bearing moment mu_max, mu_min and static system
        alpha_m = self.system.alpha_m
        qs_class_erf = self.system.qs_cl_erf  # z.B. [0, 2]
        qs_class_vorh = [self.section.qs_class_n, self.section.qs_class_p]

        if min(alpha_m) == 0:
            if qs_class_vorh[1] <= qs_class_erf[1]:
                qu = self.section.mu_max/(max(alpha_m)*self.system.l_tot ** 2)
            else:
                if self.section.section_type == "rc_rec":
                    # smooth change to 0 load bearing capacity (goal: enable more efficient optimization)
                    epsilon = 1.0e-2
                    if qs_class_vorh[1] == 1:
                        shift = 0.35
                    else:
                        shift = 0.5
                    x_d = self.section.x_p/self.section.d
                    factor = 1-0.5*(1+2/np.pi*np.arctan((x_d-shift)/epsilon))
                    qu = factor * self.section.mu_max/(max(alpha_m)*self.system.l_tot ** 2)
                else:
                    qu = 0
        else:
            if qs_class_vorh[0] <= qs_class_erf[0] & qs_class_vorh[1] <= qs_class_erf[1]:
                qu = min(self.section.mu_max/(max(alpha_m)*self.system.l_tot ** 2), self.section.mu_min /
                         (min(alpha_m)*self.system.l_tot ** 2))
            else:
                qu = 0
        return qu

    def calc_qk_zul_gzt(self, gamma_g=1.35, gamma_q=1.5):
        self.qk_zul_gzt = (self.qu - gamma_g * self.gk)/gamma_q


class Requirements:
    def __init__(self, install="ductile", lw_install=350, lw_use=350, lw_app=300):
        self.install = install
        self.lw_install = lw_install
        self.lw_use = lw_use
        self.lw_app = lw_app
