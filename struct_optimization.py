from scipy.optimize import direct
import struct_analysis
from scipy.optimize import basinhopping, Bounds  # import Minimierungsfunktion aus dem SyiPy-Paket
from scipy.optimize import minimize  # import Minimierungsfunktion aus dem SyiPy-Paket
import numpy as np

class RandomDisplacementBounds(object):
    # random displacement with bounds for basinhopping optimization
    def __init__(self, xmin, xmax, stepsize=0.1):
        self.xmin = xmin
        self.xmax = xmax
        self.stepsize = stepsize

    def __call__(self, x):
        """take a random step but ensure the new position is within the bounds """
        min_step = np.maximum(self.xmin - x, -self.stepsize)
        max_step = np.minimum(self.xmax - x, self.stepsize)

        random_step = np.random.uniform(low=min_step, high=max_step, size=x.shape)
        xnew = x + random_step

        return xnew

# OPTIMIZATION OF CROSS-SECTIONS FOR DEFINED MEMBERS
# ----------------------------------------------------------------------------------------------------------------------
def opt_rc_rec(m, to_opt="GWP", criterion="ULS", max_iter=100):
    # definition of initial values for variables, which are going to be optimized
    h0 = m.section.h  # start value for height corresponds to 1/20 of system length
    di_xu0 = m.section.bw[0][0]  # start value for rebar diameter 40 mm
    var0 = [h0, di_xu0]

    # define bounds of variables
    bh = (0.06, 1.2)  # height between 6 cm and 1.0 m
    bdi_xu = (0.006, 0.04)  # diameter of rebars between 6 mm and 40 mm
    bounds = [bh, bdi_xu]

    # definition of fixed values of cross-section
    b = m.section.b
    s_xu, di_xo, s_xo = m.section.bw[0][1], m.section.bw[1][0], m.section.bw[1][1]
    di_bw, s_bw, n_bw = m.section.bw_bg[0], m.section.bw_bg[1], m.section.bw_bg[2]
    phi, c_nom, xi, jnt_srch = m.section.phi, m.section.c_nom, m.section.xi, m.section.joint_surcharge

    co, st = m.section.concrete_type, m.section.rebar_type
    add_arg = [m.system, co, st, b, s_xu, di_xo, s_xo, m.floorstruc, m.requirements, to_opt, criterion, m.g2k, m.qk]

    # optimize with basinhopping algorithm with bounds also implemented on both levels (inner and outer):
    bounded_step = RandomDisplacementBounds(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds]))
    opt = basinhopping(rc_rqs, var0, niter=max_iter, T=1, minimizer_kwargs={"args": (add_arg,), "bounds": bounds,
                                                                            "method": "Powell"}, take_step=bounded_step)
    h, di_xu = opt.x
    optimized_section = struct_analysis.RectangularConcrete(co, st, b, h, di_xu, s_xu, di_xo, s_xo, di_bw, s_bw, n_bw,
                                                            phi, c_nom, xi, jnt_srch)
    return optimized_section

# inner function for optimizing reinforced concrete section for criteria ULS or SLS1 in terms of GWP or height
def rc_rqs(var, add_arg):
    # input: variables, which have to be optimized, additional info about cross-section and system, optimizing option
    # output: if criterion == GWP -> co2 of cross-section, punished by delta 10*(qk_zul-qk)
    # output: if criterion == h -> height of cross-section, punished by delta 1*(qk_zul-qk)
    h, di_xu = var
    system = add_arg[0]
    concrete = add_arg[1]
    reinfsteel = add_arg[2]
    b = add_arg[3]
    s_xu, di_xo, s_xo = add_arg[4:7]
    floorstruc = add_arg[7]
    criteria = add_arg[8]
    to_opt = add_arg[9]
    criterion = add_arg[10]
    g2k = add_arg[11]
    qk = add_arg[12]

    # create section
    section = struct_analysis.RectangularConcrete(concrete, reinfsteel, b, h, di_xu, s_xu, di_xo, s_xo)

    # create member
    member = struct_analysis.Member1D(section, system, floorstruc, criteria, g2k, qk)
    member.calc_qk_zul_gzt()  # calculate admissible live load

    # define penalty1, if ULS is not fulfilled
    penalty1 = max(member.qk - member.qk_zul_gzt, 0)

    # define penalty2, if SLS1 (deflections) are not fulfilled
    if member.mkd_p < member.section.mr_p and member.mkd_n < member.section.mr_n:
        d1, d2, d3 = [member.w_install - member.w_install_adm, member.w_use - member.w_use_adm,
                      member.w_app - member.w_app_adm]
    else:
        d1, d2, d3 = [member.w_install_ger - member.w_install_adm, member.w_use_ger - member.w_use_adm,
                      member.w_app_ger - member.w_app_adm]
    penalty2 = 1e5 * max(d1, d2, d3, 0)

    # define penalty3, if SLS2 (vibrations) are not fulfilled
    pen_a = member.a_ed - member.requirements.a_cd  # Grössenordnung 1e-2
    pen_w = member.wf_ed - member.requirements.w_f_cdr1 * member.r1  # HBT S. 48. r2 wird gleich 1 gesetzt
    # (Störungen im benachbarten Feld akzeptiert)  # Grössenordnung 1e-5
    pen_v = member.ve_ed - member.ve_cd  # Grössenordnung 1e-3
    if member.f1 < member.requirements.f1:
        penalty3 = max(pen_a * 1e2, pen_w * 1e5, pen_v * 1e3, 0)
    else:
        penalty3 = max(pen_w * 1e5, pen_v * 1e3, 0)

    # define penalty4, if fire resistance is not fulfilled
    member.get_fire_resistance()
    penalty4 = max(member.requirements.t_fire-member.fire_resistance, 0)

    # optimize ULS only
    if criterion == "ULS":  # optimize ultimate limit state
        if to_opt == "GWP":
            return member.section.co2*(1+penalty1)
        elif to_opt == "h":
            return member.section.h*(1+penalty1)

    # optimize SLS1 (deflections). Make sure, that also ULS is fulfilled
    elif criterion == "SLS1":  # optimize service limit state (deflections)
        if to_opt == "GWP":
            return member.section.co2*(1+penalty2)
        elif to_opt == "h":
            return member.section.h*(1+penalty2)

    # optimize SLS2 (vibrations). Make sure, that also ULS is fulfilled
    elif criterion == "SLS2":
        if to_opt == "GWP":
            to_minimize = member.section.co2*(1+penalty3)
        elif to_opt == "h":
            to_minimize = member.section.h*(1+penalty3)

    # optimize fire resistance only
    elif criterion == "FIRE":
        if to_opt == "GWP":
            return member.section.co2*(1+penalty4)
        elif to_opt == "h":
            return member.section.h * (1+penalty4)

    # optimize solution, which fulfills all requirements (ULS, SLS1 and SLS2, FIRE)
    elif criterion == "ENV":
        if to_opt == "GWP":
            to_minimize = member.section.co2*(1+penalty1+penalty2+penalty3+penalty4)
        elif to_opt == "h":
            to_minimize = member.section.h*(1+penalty1+penalty2+penalty3+penalty4)
    else:
        to_minimize = 99
        print("criterion " + criterion + " is not defined")
        print("criterion has to be 'ULS', 'SLS1', 'SLS2', 'FIRE' or 'ENV'")
    return to_minimize

# outer function for finding optimal wooden rectangular cross-section
def opt_gzt_wd_rqs(member, criterion="ULS"):
    h_0 = member.section.h
    bnds = [(0.04, 1.2)]
    minimal_h = minimize(wd_rqs_h, h_0, args=[member, criterion], bounds=bnds, method='Powell')
    h_opt = minimal_h.x[0]
    section = struct_analysis.RectangularWood(member.section.wood_type, member.section.b, h_opt)
    return section

# inner function used for optimizing wooden section in terms of height (equals co2)
def wd_rqs_h(h, args):
    m, criterion = args
    querschnitt = struct_analysis.RectangularWood(m.section.wood_type, m.section.b, h, m.section.phi)
    member = struct_analysis.Member1D(querschnitt, m.system, m.floorstruc, m.requirements, m.g2k, m.qk)
    member.calc_qk_zul_gzt()
    if criterion == "ULS":
        to_minimize = abs(member.qk - member.qk_zul_gzt)
    elif criterion == "SLS1":
        d1, d2, d3 = [member.w_install - member.w_install_adm, member.w_use - member.w_use_adm,
                      member.w_app - member.w_app_adm]
        # return penalty if w_adm =! w
        penalty2 = 1e5*max(d1, d2, d3, 0)
        to_minimize = member.section.h*(1000+penalty2)
    elif criterion == "SLS2":
        pen_a = member.a_ed - member.requirements.a_cd  # Grössenordnung 1e-2
        pen_w = member.wf_ed - member.requirements.w_f_cdr1*member.r1  # HBT S. 48. r2 wird gleich 1 gesetzt
        # (Störungen im benachbarten Feld akzeptiert)  # Grössenordnung 1e-5
        pen_v = member.ve_ed - member.ve_cd  # Grössenordnung 1e-3
        if member.f1 < member.requirements.f1:
            penalty2 = max(pen_a*1e2, pen_w*1e5, pen_v*1e3, 0)
        else:
            penalty2 = max(pen_w*1e5, pen_v*1e3, 0)
        to_minimize = member.section.h*(1+penalty2)
    elif criterion == "FIRE":
        # define penalty4, if fire resistance is not fulfilled
        member.get_fire_resistance()
        penalty4 = max(member.requirements.t_fire - member.fire_resistance, 0)
        to_minimize = member.section.h*(1+penalty4)
    elif criterion == "ENV":
        d1, d2, d3 = [member.w_install - member.w_install_adm, member.w_use - member.w_use_adm,
                      member.w_app - member.w_app_adm]
        pen_a = member.a_ed - member.requirements.a_cd  # Grössenordnung 1e-2
        pen_w = member.wf_ed - member.requirements.w_f_cdr1 * member.r1  # HBT S. 48. r2 wird gleich 1 gesetzt
        # (Störungen im benachbarten Feld akzeptiert)  # Grössenordnung 1e-5
        pen_v = member.ve_ed - member.ve_cd  # Grössenordnung 1e-3
        penalty1 = max(member.qk - member.qk_zul_gzt, 0)
        penalty2 = 1e5 * max(d1, d2, d3, 0)
        if member.f1 < member.requirements.f1:
            penalty3 = max(pen_a * 1e2, pen_w * 1e5, pen_v * 1e3, 0)
        else:
            penalty3 = max(pen_w * 1e5, pen_v * 1e3, 0)
        member.get_fire_resistance()
        penalty4 = max(member.requirements.t_fire - member.fire_resistance, 0)
        to_minimize = member.section.h * (1 + penalty1 + penalty2 + penalty3 + penalty4)
    else:
        to_minimize = 99
        print("criterion " + criterion + " is not defined")
        print("criterion has to be 'ULS', 'SLS1', 'SLS2' or ENV")
    return to_minimize


# function for returning optimal section for defined QS-type, system, requirements, loads, criterion and type of optimum
def get_optimized_section(member, criterion, to_opt, max_iter):
    if member.section.section_type == "rc_rec":
        # available to_opt arguments: "GWP", "h"
        # available criterion arguments: "ULS", "SLS1", "SLS2"
        return opt_rc_rec(member, to_opt, criterion, max_iter)
    elif member.section.section_type == "wd_rec":
        # available criterion arguments: "ULS", "SLS1", "SLS2"
        return opt_gzt_wd_rqs(member, criterion=criterion)
    elif member.section.section_type == "rc_rib":               #Added lua: does not work yet
        # available to_opt arguments: "GWP", "h"
        # available criterion arguments: "ULS", "SLS1", "SLS2"
        #return opt_rc_rib(member, to_opt, criterion, max_iter)
        print("RC-Rib optimization is being worked on")
    else:
        print("There is no optimization for the section type " + member.section.section_type + " available!")
        return member.section

# OPTIMIZATION OF CROSS-SECTIONS REGARDING PERFORMANCE WITHIN A DEFINED GWP BUDGET
# ----------------------------------------------------------------------------------------------------------------------
def get_opt_sec(section, gwp_budget):
    if section.section_type == "wd_rec":
        # outer function for finding optimal wooden rectangular cross-section
        h_0 = section.h
        bnds = [(0.02, 10.0)]
        minimal_h = minimize(wd_rec_crsc, h_0, args=[section, gwp_budget], bounds=bnds, method='Powell')
        h_opt = minimal_h.x[0]
        opt_section = struct_analysis.RectangularWood(section.wood_type, section.b, h_opt)
        return opt_section

    elif section.section_type == "rc_rec":
        # get initial values
        h_0 = section.h
        di_xu0 = section.bw[0][0]
        var0 = [h_0, di_xu0]

        # define bounds of variables
        bh = (0.06, 2.0)  # height between 6 cm and 2.0 m
        bdi_xu = (0.006, 0.04)  # diameter of rebars between 6 mm and 40 mm
        bounds = [bh, bdi_xu]

        # definition of fixed values of cross-section
        b = section.b
        s_xu, di_xo, s_xo = section.bw[0][1], section.bw[1][0], section.bw[1][1]
        di_bw, s_bw, n_bw = section.bw_bg[0], section.bw_bg[1], section.bw_bg[2]
        phi, c_nom, xi, jnt_srch = section.phi, section.c_nom, section.xi, section.joint_surcharge
        co, st = section.concrete_type, section.rebar_type
        add_arg = [co, st, b, s_xu, di_xo, s_xo, di_bw, s_bw, n_bw, phi, c_nom, xi, jnt_srch, gwp_budget]

        # optimize with basinhopping algorithm with bounds also implemented on both levels (inner and outer):
        bounded_step = RandomDisplacementBounds(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds]))
        opt = basinhopping(rc_rec_crsc, var0, minimizer_kwargs={"args": (add_arg,), "bounds": bounds,
                                                                "method": "Powell"},
                           take_step=bounded_step)
        h, di_xu = opt.x
        opt_section = struct_analysis.RectangularConcrete(co, st, b, h, di_xu, s_xu, di_xo, s_xo, di_bw, s_bw,
                                                                n_bw, phi, c_nom, xi, jnt_srch)

## XXXXXXXXXXX neuen Querschnittstyp für optimierung vorbereiten. Für mehrere parameter: basinhopping methode.

        return opt_section
    else:
        print("no optimization for section type " + section.section_type + " is defined yet within method get_opt_sec")
        return section


# inner function used for optimizing rectangular wooden section in terms of maximal bending moment and within gwp_budget
def wd_rec_crsc(h, args):
    s, gwp_budget = args
    section_updated = struct_analysis.RectangularWood(s.wood_type, s.b, h)
    penalty = 1e6*max(section_updated.co2-gwp_budget, 0)
    to_minimize = penalty - section_updated.mu_max
    return to_minimize

# inner function for optimizing reinforced concrete section in terms of maximal bending moment and within gwp_budget
def rc_rec_crsc(var, add_arg):
    h, di_xu = var
    concrete = add_arg[0]
    reinfsteel = add_arg[1]
    b = add_arg[2]
    s_xu, di_xo, s_xo = add_arg[3:6]
    di_bw, s_bw, n_bw = add_arg[6:9]
    phi, c_nom, xi, jnt_srch = add_arg[9:13]
    gwp_budget = add_arg[13]
    section_updated = struct_analysis.RectangularConcrete(concrete, reinfsteel, b, h, di_xu, s_xu, di_xo, s_xo, di_bw,
                                                          s_bw, n_bw, phi, c_nom, xi, jnt_srch)
    penalty = 1e6*max(section_updated.co2-gwp_budget, 0)
    to_minimize = penalty - section_updated.mu_max
    return to_minimize

## XXXXXXXXXXX neue funktion (returnwert wird minimiert)

