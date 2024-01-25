#!/usr/bin/env python3
import numpy as np
from scipy.interpolate import interp1d
import lx.parser

##SOME CONSTANTS##############################################
EPSILON_0 = lx.parser.EPSILON_0  # F/m
HBAR_EV = lx.parser.HBAR_EV  # eV s
HBAR_J = lx.parser.HBAR_J  # J s
MASS_E = lx.parser.MASS_E  # kg
LIGHT_SPEED = lx.parser.LIGHT_SPEED  # m/s
E_CHARGE = lx.parser.E_CHARGE  # C
BOLTZ_EV = lx.parser.BOLTZ_EV  # eV/K
AMU = lx.parser.AMU  # kg
###############################################################


##CALCULATES FLUORESCENCE LIFETIME IN S########################
def calc_lifetime(xd, yd, dyd):
    # Integrates the emission spectrum
    IntEmi = np.trapz(yd, xd)
    taxa = (1 / HBAR_EV) * IntEmi
    error = (1 / HBAR_EV) * np.sqrt(np.trapz((dyd**2), xd))
    dlife = (1 / taxa) * (error / taxa)
    return 1 / taxa, dlife


###############################################################


##CALCULATES FORSTER RADIUS####################################
def radius(xa, ya, dya, xd, yd, dyd, kappa):
    # Finds the edges of interpolation
    minA = min(xa)
    minD = min(xd)
    maxA = max(xa)
    maxD = max(xd)
    MIN = max(minA, minD)
    MAX = min(maxA, maxD)

    if MIN > MAX:
        return 0, 0
    X = np.linspace(MIN, MAX, 1000)
    f1 = interp1d(xa, ya, kind="cubic")
    f2 = interp1d(xd, yd, kind="cubic")
    f3 = interp1d(xa, dya, kind="cubic")
    f4 = interp1d(xd, dyd, kind="cubic")

    YA = f1(X)
    YD = f2(X)
    DYA = f3(X)
    DYD = f4(X)

    # Calculates the overlap
    Overlap = YA * YD / (X**4)

    # Overlap error
    OverError = Overlap * np.sqrt((DYA / YA) ** 2 + (DYD / YD) ** 2)

    # Integrates overlap
    IntOver = np.trapz(Overlap, X)

    # Integrated Overlap Error
    DeltaOver = np.sqrt(np.trapz((OverError**2), X))

    # Gets lifetime
    tau, delta_tau = calc_lifetime(xd, yd, dyd)

    # Calculates radius sixth power
    c = LIGHT_SPEED * 1e10
    const = (HBAR_EV**3) * (9 * (c**4) * (kappa**2) * tau) / (8 * np.pi)
    radius6 = const * IntOver

    # Relative error in radius6
    delta = np.sqrt((DeltaOver / IntOver) ** 2 + (delta_tau / tau) ** 2)

    # Calculates radius
    forster_radius = radius6 ** (1 / 6)

    # Error in radius
    delta_radius = forster_radius * delta / 6

    return forster_radius, delta_radius


###############################################################


##AVG HOPPING DISTANCES########################################
def r_avg(alpha, moment, rmin, dim):
    neighbors = 100
    n1 = neighbors
    n2 = 0
    n3 = 0
    if dim > 1:
        n2 = neighbors
    if dim > 2:
        n3 = neighbors
    numerator = []
    denominator = []
    for i in range(0, n1 + 1):
        for j in range(0, n2 + 1):
            for k in range(0, n3 + 1):
                r = rmin * np.sqrt(i**2 + j**2 + k**2)
                if r != 0:
                    numerator.append(r / ((alpha * moment + r) ** 6))
                    denominator.append(1.0 / ((alpha * moment + r) ** 6))
    average = np.sum(numerator) / np.sum(denominator)
    return average


###############################################################


##1D PROJECTIONS OF DIFFUSION LENGTH###########################
def difflen(forster_radius, alpha, moment, r, dim, phi):
    ld = r * (forster_radius**3) / ((alpha * moment + r) ** 3)
    ld = np.sqrt(phi) * ld / np.sqrt(dim)
    # Conversion to nm
    ld = ld / 10
    return ld


######## ANNIHILATION COEFs ###################################
def KSSA(RF, r, tau, error_radius, error_life):
    KSSA_cte = (4 * np.pi / tau) * (
        (RF) ** 6 / (r**3)
    )  # Singlet-Singlet annihilation coef
    RF_6_power_error = np.sqrt(6 * ((error_radius / RF) ** 2))
    delta = np.sqrt(RF_6_power_error**2 + (error_life / tau) ** 2)
    KSSA_error = KSSA_cte * delta
    un = 1e-24  # unit factor 1 AA^3 = 1E-24 cm^3
    KSSA_cte, KSSA_error = KSSA_cte * un, KSSA_error * un
    return KSSA_cte, KSSA_error


def KTTA(RF, r, tau, error_radius, error_life):
    KTTA_cte = (4 * np.pi / tau) * ((RF) ** 6)  # Triplet-Triplet annihilation coef
    RF_6_power_error = np.sqrt(6 * ((error_radius / RF) ** 2))
    delta = np.sqrt(RF_6_power_error**2 + (error_life / tau) ** 2)
    KTTA_error = KTTA_cte * delta
    un = 1e-48  # unit factor 1 AA^6 = 1E-48 cm^6
    KTTA_cte, KTTA_error = KTTA_cte * un, KTTA_error * un
    return KTTA_cte, KTTA_error


###############################################################


##RETURNS ABS AND EMISSION TYPES###############################
def emi_abs_types(Abs, Emi):
    abs_type = "S0"
    with open(Emi, "r",encoding='utf-8') as f:
        line = f.readlines()[1].split("->")
        emi_init = line[0].split()[-1].strip()
        emi_final = line[1].split(":")[0].strip()
    return abs_type, emi_init, emi_final


###############################################################


##Calculates average hopping distances and diffusion length for each case ####
def get_lds_dists(alpha, moment, rmin, Phi, mean_radius, error_radius):
    # convention: the last list element (with dimension 3) refers to the Amorphous one
    # dists' structure: [average distance,dimension,morphology label]
    # lds'   structure: [ld , error ld]
    dists = [
        [r_avg(alpha, moment, rmin, i), i, "%sD Crystal" % (i)] for i in range(1, 4)
    ]  # Calculates average hopping distances for each case
    dists.append([0.25 * alpha * moment + 1.25 * rmin, 3, "Amorphous "])
    lds = [
        [
            difflen(mean_radius, alpha, moment, dists[i][0], dists[i][1], Phi),
            3
            * difflen(mean_radius, alpha, moment, dists[i][0], dists[i][1], Phi)
            * (error_radius / mean_radius),
        ]
        for i in range(len(dists))
    ]
    return lds, dists


###############################################################


def run_ld(Abs, Emi, alpha, rmin, kappa, Phi):
    # Gets the Transition Dipole Moment in atomic units
    try:
        with open(Emi, "r",encoding='utf-8') as f:
            for line in f:
                moment = float(line.split("=")[1].split()[0])
                break
    except (FileNotFoundError,IndexError, ValueError):
        moment = 0

    # Loads the data
    data_abs = np.loadtxt(Abs)
    data_emi = np.loadtxt(Emi)

    # getting the types of abs and emi spectras
    abs_type, emi_init, emi_final = emi_abs_types(Abs, Emi)

    # template : abs = S0, emi = Sn ---> S0 ===> ('S0','S','S0')
    rates_types = {
        ("S0", "S", "S0"): "SSA",
        ("T1", "T", "S0"): "TTA",
        ("S0", "T", "S0"): "Only radius",
        ("S1", "S", "S0"): "Only radius",
    }

    try:
        setting_spectras = rates_types[(abs_type, emi_init[0], emi_final)]
    except KeyError:
        print(
            f"The configuration:\nabs {abs_type}\nemi {emi_init} -> {emi_final}\nis not defined! Generating standard calculations anyway ..."
        )
        setting_spectras = "Only radius"

    # Gets energies, intensities and errors for the donor and acceptor
    xa, ya, dya = data_abs[:, 0], data_abs[:, 1], data_abs[:, 2]
    xd, yd, dyd = data_emi[:, 0], data_emi[:, 1], data_emi[:, 2]

    # Lifetime calculations
    mean_life, error_life = calc_lifetime(xd, yd, dyd)

    # Radius calculations
    mean_radius, error_radius = radius(xa, ya, dya, xd, yd, dyd, kappa)

    # Diffusion length and distance calculations for all conformations
    lds, dists = get_lds_dists(alpha, moment, rmin, Phi, mean_radius, error_radius)

    with open("ld.lx", "w",encoding='utf-8') as f:
        f.write(f"Abs: {abs_type}, Emi: {emi_init} -> {emi_final}\n")
        f.write(
            f"Forster Radius:      {mean_radius:.1f} +/- {error_radius:.1f} AA \n"
        )
        f.write(
            f"Radiative Lifetime:  {mean_life:.3e} +/- {error_life:.3e} s\n"
        )
        f.write(f"Avg. Dipole Moment:  {moment:.1f} a.u. \n")

        if setting_spectras == "SSA":
            KSSA_cte, KSSA_error = KSSA(
                mean_radius, rmin, mean_life, error_radius, error_life
            )
            f.write(
                f"Annihilation Coefficient Singlet-Singlet: {KSSA_cte:5.2e} +/- {KSSA_error:5.2e} cm^3 s^-1\n"
            )
            f.write("Morphology   Avg_Hop_Distance(AA)  Diffusion_Length(nm)\n")
            for i, dist in enumerate(dists):
                f.write(f"{dist[2]}   {dist[0]:<19.1f}  {lds[i][0]:>.1f} +/- {lds[i][1]:>.1f}\n")
        if setting_spectras == "TTA":
            KTTA_cte, KTTA_error = KTTA(
                mean_radius, rmin, mean_life, error_radius, error_life
            )
            f.write(
                f"Annihilation Coefficient Triplet-Triplet: {KTTA_cte:5.2e} +/- {KTTA_error:5.2e} cm^6 s^-1\n"
            )
        if setting_spectras == "Only radius":
            pass
