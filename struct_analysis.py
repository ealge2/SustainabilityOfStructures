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
from scipy.optimize import minimize

#DEFINITONS OF MATERIAL PROPERTIES--------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class Wood:
    # defines properties of wooden material
    def __init__(self, mech_prop, database):  # retrieve basic mechanical data from database
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = ("SELECT strength_bend, strength_shea, E_modulus, density_load, burn_rate FROM material_prop WHERE"
                   " name="+mech_prop)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fmk, self.fvd, self.Emmean, self.weight, self.burn_rate = result[0]
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
    def __init__(self, mech_prop, database, dmax=32):  # retrieve basic mechanical data from database (self, table,
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
        self.dmax = dmax

    def get_design_values(self, gamma_c=1.5, eta_t=1):  # calculate design values
        eta_fc = min((30e6/self.fck) ** (1/3), 1)  # SIA 262, 4.2.1.2, Formel (26)
        self.fcd = self.fck * eta_fc * eta_t / gamma_c  # SIA 262, 2.3.2.3, Formel (2)
        self.tcd = 0.3 * eta_t * 1e6 * (self.fck*1e-6) ** 0.5/gamma_c  # SIA 262, 2.3.2.4, Formel (3)
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

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class Section:
    # contains fundamental section properties like section type weight, resistance and stiffness
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

class RectangularWood(SupStrucRectangular, Section):
    # defines properties of rectangular, wooden cross-section
    def __init__(self, wood_type, b, h, phi=0.6, xi=0.01, ei_b=0.0):  # create a rectangular timber object
        section_type = "wd_rec"
        super().__init__(section_type, b, h, phi)
        self.wood_type = wood_type
        mu_el, vu_el = self.calc_strength_elast(wood_type.fmd, wood_type.fvd)
        self.mu_max, self.mu_min = [mu_el, -mu_el]   #Readme: Why is this needed for wood? -> is not needed for wood.
        # However, as the same resistance values should be provided for all cross-sections, I defined them for both
        # directions for wood too
        self.vu_p, self.vu_n = vu_el, vu_el
        self.qs_class_n, self.qs_class_p = [3, 3]   # Required cross-section class: 1:=PP, 2:EP, 3:EE
        self.g0k = self.calc_weight(wood_type.weight)
        self.ei1 = self.wood_type.Emmean*self.iy  # elastic stiffness [Nm^2]
        self.co2 = self.a_brutt * self.wood_type.GWP * self.wood_type.density  # [kg_CO2_eq/m]
        self.cost = self.a_brutt * self.wood_type.cost
        self.ei_b = ei_b  # stiffness perpendicular to direction of span
        self.xi = xi  # damping factor, preset value see: HBT, Page 47 (higher value for some buildups possible)

    @staticmethod
    def fire_resistance(member):
        bnds = [(0, 240)]
        t0 = 60
        max_t = minimize(RectangularWood.fire_minimizer, t0, args=[member], bounds=bnds)
        t_max = max_t.x[0]
        return t_max

    @staticmethod
    def fire_minimizer(t, args):
        member = args[0]
        rem_sec = RectangularWood.remaining_section(member.section, member.fire, t)
        mu_fire = 1.8*rem_sec.mu_max
        vu_fire = 1.8*rem_sec.vu_p  # SIA 265 (51)
        qd_fire = member.psi[2]*member.qk + member.gk
        qd_fire_zul = min(mu_fire/(max(member.system.alpha_m) * member.system.l_tot**2),
                          vu_fire/(max(member.system.alpha_v) * member.system.l_tot))
        to_opt = abs(qd_fire-qd_fire_zul)
        return to_opt


    @staticmethod
    def remaining_section(section, fire, t=60, dred=0.007):
        betan = section.wood_type.burn_rate
        dcharn = betan*t
        d_ef = dcharn + dred
        h_fire = max(section.h-d_ef*(fire[0]+fire[2]))
        b_fire = max(section.b-d_ef*(fire[1]+fire[3]), 0)
        rem_sec = RectangularWood(section.wood_type, b_fire, h_fire)
        return rem_sec

class RectangularConcrete(SupStrucRectangular):
    # defines properties of rectangular, reinforced concrete cross-section
    def __init__(self, concrete_type, rebar_type, b, h, di_xu, s_xu, di_xo, s_xo, di_bw=0.0, s_bw=0.15, n_bw=0,
                 phi=2.0, c_nom=0.03, xi=0.02):

        # create a rectangular concrete object
        section_type = "rc_rec"
        super().__init__(section_type, b, h, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom
        self.bw = [[di_xu, s_xu], [di_xo, s_xo]]
        self.bw_bg = [di_bw, s_bw, n_bw]
        mr = self.b * self.h ** 2 / 6 * 1.3 * self.concrete_type.fctm  # cracking moment
        self.mr_p, self.mr_n = mr, mr
        [self.d, self.ds] = self.calc_d()
        [self.mu_max, self.x_p, self.as_p, self.qs_class_p] = self.calc_mu('pos')
        [self.mu_min, self.x_n, self.as_n, self.qs_class_n] = self.calc_mu('neg')
        self.roh, self.rohs = self.as_p/self.d, self.as_n/self.ds
        [self.vu_p, self.vu_n, self.as_bw] = self.calc_shear_resistance()
        self.g0k = self.calc_weight(concrete_type.weight)
        a_s_tot = self.as_p + self.as_n + self.as_bw
        co2_rebar = a_s_tot * self.rebar_type.GWP * self.rebar_type.density  # [kg_CO2_eq/m]
        co2_concrete = (self.a_brutt-a_s_tot) * self.concrete_type.GWP * self.concrete_type.density  # [kg_CO2_eq/m]
        self.ei1 = self.concrete_type.Ecm*self.iy  # elastic stiffness concrete (uncracked behaviour) [Nm^2]
        self.co2 = co2_rebar + co2_concrete
        self.cost = (a_s_tot * self.rebar_type.cost + (self.a_brutt-a_s_tot) * self.concrete_type.cost
                     + self.concrete_type.cost2)
        self.ei_b = self.ei1
        self.xi = xi  # XXXXXXX preset value is an assumption. Has to be verified with literature. XXXXXXX
        self.ei2 = self.ei1 / self.f_w_ger(self.roh, self.rohs, 0, self.h, self.d)

    def calc_d(self):
        d = self.h - self.c_nom - self.bw[0][0]/2
        ds = self.h - self.c_nom - self.bw[1][0]/2
        return d, ds

    def calc_mu(self, sign='pos'):
        b = self.b
        fsd = self.rebar_type.fsd
        fcd = self.concrete_type.fcd
        if sign == 'pos':
            [mu, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[0][0], self.bw[0][1], self.d, b, fsd, fcd, self.mr_p)
        elif sign == 'neg':
            [mus, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[1][0], self.bw[1][1], self.ds, b, fsd, fcd, self.mr_n)
            mu = -mus
        else:
            [mu, x, a_s, qs_klasse] = [0, 0, 0, 0]
            print("sigen of moment resistance has to be 'neg' or 'pos'")
        return mu, x, a_s, qs_klasse


    @staticmethod
    def mu_unsigned(di, s, d, b, fsd, fcd, mr):
        # units input: [m, m, m, m, N/m^2, N/m^2]
        a_s = np.pi * di ** 2 / (4 * s) * b  # [m^2]
        omega = a_s * fsd / (d * b * fcd)  # [-]
        mu = a_s * fsd * d * (1-omega/2)  # [Nm]
        x = omega * d / 0.85  # [m]
        if x/d <= 0.35 and mu >= mr:
            return mu, x, a_s, 1
        elif x/d <= 0.5 and mu >= mr:
            return mu, x, a_s, 2
        else:
            return mu, x, a_s, 99  # Querschnitt hat ungenügendes Verformungsvermögen


    def calc_shear_resistance(self, d_installation=0.0):
        # calculates shear resistance with d
        di = self.bw_bg[0]  # diameter
        s = self.bw_bg[1]  # spacing
        n = self.bw_bg[2]  # number of stirrups per spacing
        fck = self. concrete_type.fck
        fcd = self.concrete_type.fcd
        tcd = self.concrete_type.tcd
        dmax = self.concrete_type.dmax  # dmax in mm
        fsk = self.rebar_type.fsk
        fsd = self.rebar_type.fsd
        es = self.rebar_type.Es
        bw = self.b
        d = self.d
        ds = self.ds
        x_p = self.x_p
        x_n = self.x_n
        as_bw = np.pi*di**2/4*n/s
        if d_installation < d/6:
            dv_p = d
        else:
            dv_p = d-d_installation
        if d_installation < ds/6:
            dv_n = ds
        else:
            dv_n = ds-d_installation
        vu_p = self.vu_unsigned(bw, as_bw, d, dv_p, x_p, fck, fcd, tcd, fsk, fsd, es, dmax)
        vu_n = self.vu_unsigned(bw, as_bw, ds, dv_n, x_n, fck, fcd, tcd, fsk, fsd, es, dmax)
        return vu_p, vu_n, as_bw

    @staticmethod
    def vu_unsigned(bw, as_bw, d, dv, x, fck, fcd, tcd, fsk, fsd, es, dmax=32, alpha=np.pi/4, kc=0.55):
        if as_bw == 0:  # cross-section without stirrups
            ev = 1.5*fsd/es
            kg = 48/(16+dmax)
            kd = 1/(1+ev*d*kg)
            vrd = kd*tcd*dv
            return vrd
        else:  # cross-section with vertical stirrups
            z = d-0.85*x/2
            vrds = as_bw*z*fsd
            vrdc = bw*z*kc*fcd*np.sin(alpha)*np.cos(alpha)  # unit of alpha: [rad]
            rohw = as_bw/bw
            rohw_min = 0.001*(fck*1e-6/30)**0.5*500/(fsk*1e-6)
            if rohw < rohw_min:
                print("minimal reinforcement ratio of stirrups is lower than required according to SIA 262, (110)")
            return min(vrds, vrdc)

    @staticmethod
    def f_w_ger(roh, rohs, phi, h, d):
        f = (1-20*rohs)/(10*roh**0.7)*(0.75+0.1*phi)*(h/d)**3
        return f

    @staticmethod
    def fire_resistance(section):
        # fire resistance of 1-D load-bearing plates according to SIA 262, Tab.16
        c_nom = section.c_nom
        h = section.h
        b = section.b
        if c_nom >= 0.04 and h >= 0.15 and b >= 0.4:
            resistance = 180
        elif c_nom >= 0.03 and h >= 0.12 and b >= 0.3:
            resistance = 120
        elif c_nom >= 0.03 and h >= 0.1 and b >= 0.2:
            resistance = 90
        elif c_nom >= 0.02 and h >= 0.08 and b >= 0.15:
            resistance = 60
        elif c_nom >= 0.02 and h >= 0.06 and b >= 0.1:
            resistance = 30
        else:
            resistance = 0
        return resistance

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class SupStrucRipped(Section):

    def __init__(self,section_type,b,b_w,h,h_r,phi=0):
        super().__init__(section_type)
        self.b = b # flange width [m] (Abstand Rippenachse-Rippenachse)
        self.b_w = b_w #web width [m]


    def calc_area(self):


    def calc_beff(self):

    def calc_moment_of_inertia(self):

    def calc_strength_elast(self, fy, ty):

    def calc_strength_plast(self, fy, ty):

    def calc_weight(self,spec_weight):




class RippedConcrete(SupStrucRipped):
    # defines properties of a rectangular, reinforced concrete section
    def __init__(self, concrete_type, rebar_type, b, b_w, h, h_f, di_xu, s_xu, di_xo, s_xo, di_xw, n_xw, di_bw, s_bw, n_bw=0
                 phi=2.0, c_nom=0.03, xi=0.02):
        section_type = "rc_rec"
        super().__init__(section_type, b, b_w, h, h_f, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom
        self.bw 0 =
        self.bw_bg
        mr =

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class MatLayer:  # create a material layer
    def __init__(self, mat_name, h_input, roh_input, database):  # get initial data from database
        self.name = mat_name
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get properties from database
        inquiry = "SELECT h_fix, E, density, weight, GWP FROM floor_struc_prop WHERE name=" + mat_name
        cursor.execute(inquiry)
        result = cursor.fetchall()
        h_fix, e, density, weight, self.GWP = result[0]
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
        if e == None:
            self.ei = 0.0
        else:
            i = 1 * self.h**3 / 12
            self.ei = e*i
        self.gk = self.weight * self.h  # weight per area in N/m^2
        self.co2 = self.density * self.h * self.GWP  # CO2-eq per area in kg-C02/m^2


class FloorStruc:  # create a floor structure
    def __init__(self, mat_layers, database_name):
        self.layers = []
        self.co2 = 0
        self.gk_area = 0
        self.h = 0
        self.ei = 0
        for mat_name, h_input, roh_input in mat_layers:
            current_layer = MatLayer(mat_name, h_input, roh_input, database_name)
            self.layers.append(current_layer)
            self.co2 += current_layer.co2
            self.gk_area += current_layer.gk
            self.h += current_layer.h
            self.ei = max(self.ei, current_layer.ei)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class BeamSimpleSup:
    def __init__(self, length):
        self.l_tot = length
        self.li_max = self.l_tot  # max span (used for calculation of admissible deflections)
        self.alpha_m = [0, 1/8]  # Faktor zur Berechung des Momentes unter verteilter Last
        self.alpha_v = [0, 1/2] # Faktor zur Berechung der Querkarft unter verteilter Last
        self.qs_cl_erf = [3, 3]  # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 5/384  # Faktor zur Berechung der Durchbiegung unter verteilter Last
        self.kf2 = 1.0  # Hilfsfaktor zur Brücksichtigung der Spannweitenverhältnisse bei Berechnung f1 gem. HBT, S. 46
        self.alpha_w_f_cd = 1/48  # Faktor zur Berechung der Durchbiegung unter Einzellast


class Member1D:
    def __init__(self, section, system, floorstruc, requirements, g2k=0.0, qk=2e3, psi0=0.7, psi1=0.5, psi2=0.3,
                 fire_b=True, fire_l=False, fire_t=False, fire_r=False):
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
        self.m = self.q_per/10
        self.w_install_adm = self.system.li_max/self.requirements.lw_install
        self.w_use_adm = self.system.li_max/self.requirements.lw_use
        self.w_app_adm = self.system.li_max/self.requirements.lw_app
        self.qu = self.calc_qu()
        self.mkd_n = self.system.alpha_m[0] * (self.gk + self.qk) * self.system.l_tot**2
        self.mkd_p = self.system.alpha_m[1] * (self.gk + self.qk) * self.system.l_tot ** 2
        self.qk_zul_gzt = float
        self.fire = [0, 0, 0, 0]  # fire from bottom, left, top, right (0: no fire; 1: fire)
        if fire_b is True:
            self.fire[0] = 1
        if fire_l is True:
            self.fire[1] = 1
        if fire_t is True:
            self.fire[2] = 1
        if fire_r is True:
            self.fire[3] = 1
        self.fire_resistance = []

        # calculation of deflections (uncracked plus cracked for concrete sections)
        section_material = self.section.section_type[0:2]
        unit_def = self.system.alpha_w * self.system.l_tot ** 4 / self.section.ei1  # deflection for q = 1, phi = 0
        if self.requirements.install == "ductile":
            self.w_install = unit_def * (self.q_freq + self.q_per * (self.section.phi - 1))
            if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                self.w_install_ger = unit_def * (
                    self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi,
                                                             self.section.h, self.section.d)
                    + (self.q_freq - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,
                                                                               0, self.section.h, self.section.d)
                    - self.q_per
                )
        elif self.requirements.install == "brittle":
            self.w_install = unit_def * (self.q_rare + self.q_per * (self.section.phi - 1))
            if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                self.w_install_ger = unit_def * (
                    self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi,
                                                             self.section.h, self.section.d)
                    + (self.q_rare - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,
                                                                               0, self.section.h, self.section.d)
                    - self.q_per
                )
        self.w_use = unit_def * (self.q_freq - self.gk)
        if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
            self.w_use_ger = unit_def * (
                (self.q_freq - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, 0,
                                                                         self.section.h, self.section.d)
            )
        self.w_app = unit_def * (self.q_per * (1 + self.section.phi))
        if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
            self.w_app_ger = unit_def * (
                    self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi,
                                                             self.section.h, self.section.d)
            )
        self.co2 = system.l_tot * (self.floorstruc.co2 + self.section.co2)

        # calculation first frequency (uncracked cross-section, method for cracked cross-section is not implemented jet)
        self.f1 = self.calc_f1()
        # calculation of further vibration criteria for wooden cross-sections
        section_material = self.section.section_type[0:2]
        if section_material == "wd" or section_material == "rc":  # check for material type
            self.ei_b = max(self.section.ei_b,
                            self.floorstruc.ei)  # Berücksichtigung n.t. Bodenaufbau gemäss Beispielsammlung HBT)
            self.bm_rech = self.system.li_max / 1.1 * (self.ei_b / self.section.ei1) ** 0.25  # HBT Seite 46
            self.a_ed = self.calc_vib1()
            self.wf_ed, self.ve_ed = self.calc_vib2()
            if self.section.xi < 0.015:
                self.r1 = 1.0  # HBT S. 48
            elif self.section.xi < 0.025:
                self.r1 = 1.15  # HBT S. 48
            else:
                self.r1 = 1.25  # HBT S. 48
            self.ve_cd = self.requirements.alpha_ve_cd*100**(self.f1*self.section.xi-1)

    def calc_qu(self):
        # calculates maximal load qu in respect to bearing moment mu_max, mu_min and static system
        alpha_m = self.system.alpha_m
        alpha_v = self.system.alpha_v
        qs_class_erf = self.system.qs_cl_erf  # z.B. [0, 2]
        qs_class_vorh = [self.section.qs_class_n, self.section.qs_class_p]

        if min(alpha_m) == 0:
            if qs_class_vorh[1] <= qs_class_erf[1]:
                qu_bend = self.section.mu_max/(max(alpha_m)*self.system.l_tot ** 2)
            else:
                if self.section.section_type == "rc_rec":
                    # smooth change to 0 load bearing capacity when roh<roh_min or roh>roh_zul
                    # (enables more efficient optimization)
                    epsilon = 1.0e-3
                    if qs_class_vorh[1] == 1:
                        shift = 0.35
                    else:
                        shift = 0.5
                    x_d = self.section.x_p/self.section.d
                    factor = min(0.5 * (1 + 2 / np.pi * np.arctan((self.section.mu_max-self.section.mr_p) / epsilon)),
                                        1-0.5*(1+2/np.pi*np.arctan((x_d-shift)/epsilon)))
                    qu_bend = factor * self.section.mu_max/(max(alpha_m)*self.system.l_tot ** 2)
                else:
                    qu_bend = 0
            qu_shear = self.section.vu_p / (max(alpha_v) * self.system.l_tot)
        else:
            if qs_class_vorh[0] <= qs_class_erf[0] & qs_class_vorh[1] <= qs_class_erf[1]:
                qu_bend = min(self.section.mu_max/(max(alpha_m)*self.system.l_tot ** 2), self.section.mu_min /
                              (min(alpha_m)*self.system.l_tot ** 2))
            else:
                qu_bend = 0
            qu_shear = min(self.section.vu_p/(max(alpha_v) * self.system.l_tot),
                           self.section.vu_n/(min(alpha_v)*self.system.l_tot))
        return min(qu_bend, qu_shear)

    def calc_qk_zul_gzt(self, gamma_g=1.35, gamma_q=1.5):
        self.qk_zul_gzt = (self.qu - gamma_g * self.gk)/gamma_q

    def calc_f1(self):
        # calculates first frequency of system according to HBT, Seite 46
        kf2 = self.system.kf2
        l_rech = self.system.li_max
        section_material = self.section.section_type[0:2]
        if section_material == "rc":  # take cracked stiffness for calculation of concrete sections if section is cracked
            if self.mkd_p < self.section.mr_p and self.mkd_n > self.section.mr_n:
                eil = self.section.ei1
            else:
                eil = self.section.ei2
        else:
            eil = self.section.ei1
        m = self.m
        f1 = kf2*np.pi/(2*l_rech**2)*(eil/m)**0.5  # HBT, Seite 46
        return f1

    def calc_vib1(self, f0=700):
        # calculates a_Ed according to HBT, Seite 47
        f1 = self.f1
        m_gen = self.m * self.system.li_max / 2 * self.bm_rech
        xi = self.section.xi
        if f1 <= 5.1:
            alpha = 0.2
            ff = f1
        elif f1 <= 6.9:
            alpha = 0.06
            ff = f1
        else:
            alpha = 0.06
            ff = 6.9
        a_ed = 0.4*f0*alpha/m_gen*1/(((f1/ff)**2-1)**2+(2*xi*f1/ff)**2)**0.5  # HBT, Seite 47
        return a_ed

    def calc_vib2(self, f=1000):
        # calculates W_F,ED according to to HBT, Seite 48
        wf_ed = self.system.alpha_w_f_cd*f*self.system.li_max**3/(self.bm_rech*self.section.ei1)
        section_material = self.section.section_type[0:2]
        if section_material == "rc":  # take cracked stiffness for calculation of concrete sections
            eil = self.section.ei2
        else:
            eil = self.section.ei1
        ve_ed = 364/(self.bm_rech*(self.m**3*eil*1e6)**0.25)
        return wf_ed, ve_ed

    def get_fire_resistance(self):
        # evaluate fire resistance
        if self.section.section_type == "rc_rec":
            fire_resistance = RectangularConcrete.fire_resistance(self.section)
        elif self.section.section_type == "wd_rec":
            fire_resistance = RectangularWood.fire_resistance(self)
        else:
            print("fire resistance for is not defined for that cross-section type.")
            fire_resistance = None
        self.fire_resistance = fire_resistance
class Requirements:
    def __init__(self, install="ductile", lw_install=350, lw_use=350, lw_app=300, f1=8, a_cd=0.1, w_f_cdr1=1.0e-3,
                 alpha_ve_cd=1/3, fire='R60'):
        self.install = install
        self.lw_install = lw_install  # preset value: SIA 260
        self.lw_use = lw_use  # preset value: SIA 260
        self.lw_app = lw_app  # preset value: SIA 260
        self.f1 = f1  # preset value: HBT, Seite 46
        self.a_cd = a_cd  # preset value: HBT, Seite 46
        self.w_f_cdr1 = w_f_cdr1  # preset value: HBT, Seite 48
        self.alpha_ve_cd = alpha_ve_cd  # preset value: HBT, Seite 49
        self.t_fire = int(fire[1:])  # unit: [min]
