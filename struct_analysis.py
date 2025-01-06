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

#DEFINITONS OF MATERIAL PROPERTIES--------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
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
        self.mu_max, self.mu_min = [mu_el, mu_el]   #Readme: Why is this needed for wood?
        self.vu = vu_el
        self.qs_class_n, self.qs_class_p = [3, 3]   #Readme: What is this needed for?
        self.g0k = self.calc_weight(wood_type.weight)
        self.ei1 = self.wood_type.Emmean*self.iy  # elastic stiffness [Nm^2]
        self.co2 = self.a_brutt * self.wood_type.GWP * self.wood_type.density  # [kg_CO2_eq/m]
        self.cost = self.a_brutt * self.wood_type.cost
        self.xi = xi  # damping factor, preset value see: HBT, Page 47 (higher value for some buildups possible)
        self.ei_b = ei_b  # stiffness perpendicular to direction of span


class RectangularConcrete(SupStrucRectangular):
    # defines properties of rectangular, reinforced concrete cross-section
    def __init__(self, concrete_type, rebar_type, b, h, di_xu, s_xu, di_xo, s_xo, di_bg, s_bg_l, s_bg_t, phi=2.0, c_nom=0.03):

        # create a rectangular concrete object
        section_type = "rc_rec"
        super().__init__(section_type, b, h, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom
        self.bw = [[di_xu, s_xu], [di_xo, s_xo]]
        self.bw_bg = [di_bg, s_bg_l, s_bg_t]    #diameter, distance longitudinal, distance transversal
        [self.d, self.ds] = self.calc_d()
        [self.mu_max, self.x_p, self.as_p, self.qs_class_p] = self.calc_mu('pos')
        [self.mu_min, self.x_n, self.as_n, self.s_class_n] = self.calc_mu('neg')
        [self.vu, self.as_bg] = self.calc_vu()
        self.g0k = self.calc_weight(concrete_type.weight)
        a_s_tot = self.as_p + self.as_n  + self.as_bg
        co2_rebar = a_s_tot * self.rebar_type.GWP * self.rebar_type.density  # [kg_CO2_eq/m]
        co2_concrete = (self.a_brutt-a_s_tot) * self.concrete_type.GWP * self.concrete_type.density  # [kg_CO2_eq/m]
        self.ei1 = self.concrete_type.Ecm*self.iy  # elastic stiffness concrete (uncracked behaviour) [Nm^2]
        self.co2 = co2_rebar + co2_concrete
        self.cost = (a_s_tot * self.rebar_type.cost + (self.a_brutt-a_s_tot) * self.concrete_type.cost
                     + self.concrete_type.cost2)
        # self.xi = XX  no damping factor for concrete defined (no value needed for applied calculation concepts)
        # self.ei_b = XX  no stiffness in perpendicular direction defined (no value needed for applied calculation concepts)
        #   self.ei2 = # XXXXXXXXXXToDoXXXXXXXXXX

    def calc_d(self):
        d = self.h - self.c_nom - self.bw[0][0]/2
        ds = self.h - self.c_nom - self.bw[1][0]/2
        return d, ds

    def calc_mu(self, sign='pos'):
        b = self.b, fsd = self.rebar_type.fsd, fcd = self.concrete_type.fcd
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

    def calc_vu(self, D_max=32, alpha = 45):
        b = self.b
        z = 0.9 * self.d #APPROXIMATION
        as_bg = np.pi * self.bw_bg[0] ** 2 / (4 * self.bw_bg[1] * self.bw_bg[2]) * b
        if self.bw_bg[0] == 0:   #Bauteile ohne Querkraftbewehrung
            k_g = 48/(16+D_max)
            e_v = 1.5*self.rebar_type.fsd*self.rebar_type.Es #(39) -> READ ME: Überlegen, wie Formel (38) implementiert wird
            k_d = 1/(1+e_v*self.d*k_g)
            vu_c = k_d*self.concrete_type.tcd*self.d
            vu_s = 0
            vu = vu_c
            return vu_c, vu_s, vu , as_bg
        else:
            vu_s = as_bg*z*self.rebar_type.fsd*1/np.tan(np.radians(alpha)) #Einschnittige Bügel gerechnet
            k_c = 0.55 #READ ME: make kc variable?
            vu_c = b*z*k_c*self.concrete_type.fcd*np.sin(np.radians(alpha))*np.cos(np.radians(alpha))
            vu = min(vu_c,vu_s)
            return vu_c, vu_s, vu, as_bg

#Ripped cross sections
class SupStrucRipped(Section):
    # defines cross-section dimensions and has methods to calculate static properties of ribbed,
    # non-cracked sections
    def __init__(self, section_type, b, b_w, h, h_f, phi=0):  # create a rectangular object
        super().__init__(section_type)
        self.b_w = b_w  # web width [m]
        self.b = b      # width mid rib to mid rib [m]
        self.h = h      # total height [m]
        self.h_f = h_f  # flange height [m]
        #self.a_brutt = self.calc_area()
        #self.iy = self.calc_moment_of_inertia()
        self.phi = phi

        def calc_area(self):
            #  in: width b [m], web with b_w [m], total height h [m], flange height h_f [m]
            #  out: area [m^2]
            a_brutt = self.b * self.h_f + self.b_w*(self.h-self.h_f)
            return a_brutt

class RippedConcrete(SupStrucRipped):
    # defines properties of a rectangular, reinforced concrete section

    def __init__(self, concrete_type, rebar_type, b, b_w, h, h_f, di_xu, s_xu, di_xo, s_xo, di_xw, n_xw, di_bg, s_bg, l0, phi=2.0, c_nom=0.03):
        section_type = "rc_rec"
        super().__init__(section_type, b, b_w, h, h_f, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom
        self.l0 = l0
        self.bw = [[di_xw, n_xw], [di_xu, s_xu], [di_xo, s_xo]]
        self.bw_bg = [di_bg, s_bg]
        [self.d, self.dso, self.dsu] = self.calc_d()
        self.zs = self.calc_zs()
        self.b_eff = self.calc_beff()
        self.iy = self.calc_moment_of_inertia()

        [self.mu_max, self.x_p, self.as_p, self.qs_class_p] = self.calc_mu('pos')
        [self.mu_min, self.x_n, self.as_n, self.qs_class_n] = self.calc_mu('neg')

    def calc_beff(self):
        #computes effective width of concrete flange: SIA 262 4.1.3.3.2
        b = self.b, b_w = self.b_w, l0 = self.l0
        b_effi = min(0.2*(b-b_w)/2+0.1*l0, 0.2*l0)  #SIA 262 (20)
        b_eff = min(2*b_effi + b_w, b)              #SIA 262 (19)
        return b_eff

    def calc_d(self):
        #calculates
        d = self.h-self.c_nom-self.bw_bg [1][1]-self.bw[0][0]/2
        dso = self.h_f-self.c_nom -self.bw[1][1]/2
        dsu = self.h_f-self.c_nom-self.bw[2][2]/2
        return d, dso, dsu

    def calc_zs(self):
        b = self.b, b_w = self.b_w, h = self.h, h_f = self.h_f
        zs = (b_w*(h-h_f)*(h-h_f)/2+b*h_f*h_f/2)/(b_w*(h-h_f)+b*h_f)
        return zs

    def calc_moment_of_inertia(self):
        #  in: width b [m], height h [m]
        #  out: second moment of inertia Iy [m^4]
        iy_rib = self.b_w* (self.h-self.h_f) ** 3 / 12
        iy_flange = self.b* self.h_f ** 3 / 12
        sa_rib = self.b_w*(self.h-self.h_f)*abs(self.zs-(self.h-self.h_f))**2
        sa_flange = self.b*self.h*abs(self.zs-self.h/2)**2
        iy = iy_rib + iy_flange + sa_rib + sa_flange
        return iy

    def calc_mu(self, sign='pos'):
        b_w = self.b_w
        b_eff = self.b_eff
        fsd = self.rebar_type.fsd
        fcd = self.concrete_type.fcd
        if sign == 'pos':
            [mu, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[0][0], self.bw[0][1], self.d, b_eff, fsd, fcd)
        elif sign == 'neg':
            [mu, x, a_s, qs_klasse] = print("to implement") #self.mu_unsigned(self.bw[1][0], self.bw[1][1], self.dso, b_w, fsd, fcd)
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
        if x > self.h_fl
            print("x>h_fl")
        else:
            print("x<hfl")
        if x/d <= 0.35:
            return mu, x, a_s, 1
        elif x/d <= 0.5:
            return mu, x, a_s, 2
        else:
            return mu, x, a_s, 99  # Querschnitt hat ungenügendes Verformungsvermögen


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
            self.ei = 0
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
        self.alpha_m = [0, 1/8]
        self.qs_cl_erf = [3, 3]  # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 5/384
        self.kf2 = 1.0  # Hilfsfaktor zur Brücksichtigung der Spannweitenverhältnisse bei Berechnung f1 gem. HBT, S. 46
        self.alpha_w_f_cd = 1/48


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
        self.m = self.q_per/10
        self.w_install_adm = self.system.li_max/self.requirements.lw_install
        self.w_use_adm = self.system.li_max/self.requirements.lw_use
        self.w_app_adm = self.system.li_max/self.requirements.lw_app
        self.qu = self.calc_qu()
        self.qk_zul_gzt = float
        self.ei_b = max(self.section.ei_b,
                   self.floorstruc.ei)  # Berücksichtigung n.t. Bodenaufbau gemäss Beispielsammlung HBT)
        self.bm_rech = self.system.li_max / 1.1 * (self.ei_b / self.section.ei1) ** 0.25  # HBT Seite 46

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

        # calculation first frequency (uncracked cross-section, method for cracked cross-section is not implemented jet)
        self.f1 = self.calc_f1()
        # calculation of further vibration criteria for wooden cross-sections
        section_material = self.section.section_type[0:2]
        if section_material == "wd":  # check for material type
            if (self.f1 < 8.0):  # check for frequency below 8 Hz
                self.a_ed = self.calc_vib1()
            self.wf_ed, self.ve_ed = self.calc_vib2()



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

    def calc_f1(self):
        # calculates first frequency of system according to HBT, Seite 46
        kf2 = self.system.kf2
        l_rech = self.system.li_max
        eil = self.section.ei1
        m = self.m
        f1 = kf2*np.pi/(2*l_rech**2)*(eil/m)**0.5  # HBT, Seite 46
        return f1

    def calc_vib1(self, f0=700):
        # calculates a_Ed according to HBT, Seite 47
        f1 = self.f1
        if f1 >= 8:
            print("frequency f1 above 8 Hz, no evaluation of a_Ed needed")
            return 0
        else:
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
        ve_ed = 367/(self.bm_rech*(self.m**3*self.section.ei1*1e6)**0.25)
        return wf_ed, ve_ed


class Requirements:
    def __init__(self, install="ductile", lw_install=350, lw_use=350, lw_app=300, f1=8, a_cd=0.1, w_f_cdr1=1.0):
        self.install = install
        self.lw_install = lw_install  # preset value: SIA 260
        self.lw_use = lw_use  # preset value: SIA 260
        self.lw_app = lw_app  # preset value: SIA 260
        self.f1 = f1  # preset value: HBT, Seite 46
        self.a_cd = a_cd  # preset value: HBT, Seite 46
        self.w_f_cdr1 = w_f_cdr1  # preset value: HBT, Seite 48
