import sys
import os
import re
import numpy as np

ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca "
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr "
    "Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
    "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm "
    "Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def element_symbol(atom):
    """Return a chemical symbol for an element label or atomic number."""
    if atom.isdigit():
        number = int(atom)
        if not 1 <= number < len(ELEMENTS):
            raise ValueError("Unsupported atomic number: " + atom)
        return ELEMENTS[number]
    symbol = atom.capitalize()
    if symbol not in ELEMENTS[1:]:
        raise ValueError("Unsupported element label: " + atom)
    return symbol


def read_gaussian_input(filename):
    """Read a single-job Gaussian input with an explicit Cartesian geometry.

    Preserve route, Link 0 settings, title, and post-geometry sections, including
    custom basis sets and solvent data. Checkpoint-only and Z-matrix inputs are
    rejected rather than interpreted as Cartesian coordinates.
    """
    with open(filename, encoding="utf-8") as handle:
        text = handle.read()
    if "--link1--" in text.lower():
        raise ValueError("Use a single-job Gaussian template; opt/freq Link1 is generated automatically.")
    sections = re.split(r"\n\s*\n", text.strip())
    if len(sections) < 3:
        raise ValueError("Missing Gaussian route, title, or Cartesian geometry.")
    header = sections[0].splitlines()
    route_start = next((i for i, line in enumerate(header) if line.lstrip().startswith("#")), None)
    if route_start is None:
        raise ValueError("Missing Gaussian route section.")
    link0 = [line.strip() for line in header[:route_start] if line.strip()]
    if any(not line.startswith("%") for line in link0):
        raise ValueError("Unexpected text before the Gaussian route.")
    route = " ".join(line.strip() for line in header[route_start:])
    if re.search(r"\b(?:geom\s*=\s*(?:\([^)]*)?(?:allcheck|check\w*)|guess\s*=\s*(?:\([^)]*)?read)\b", route, re.I):
        raise ValueError("The starting template must not depend on a checkpoint geometry or wavefunction.")
    if any(line.lower().startswith("%oldchk") for line in link0):
        raise ValueError("%oldchk is not supported in a fresh conformer template.")
    molecule = sections[2].splitlines()
    cm = molecule[0].split()
    if len(cm) != 2:
        raise ValueError("Expected one charge and multiplicity pair.")
    charge, multiplicity = map(int, cm)
    if multiplicity < 1:
        raise ValueError("Multiplicity must be positive.")
    atoms, coordinates = [], []
    for line in molecule[1:]:
        fields = line.replace(",", " ").split()
        if len(fields) != 4:
            raise ValueError("Expected element X Y Z; frozen atoms and Z-matrices are not supported.")
        atoms.append(element_symbol(fields[0]))
        coordinates.append([float(value.replace("D", "E").replace("d", "e")) for value in fields[1:]])
    if not atoms or not np.isfinite(coordinates).all():
        raise ValueError("Missing or invalid Cartesian coordinates.")
    nproc, mem = 1, "1GB"
    for line in link0:
        name, _, value = line.partition("=")
        name = name.strip()
        if name.lower() in ("%nproc", "%nprocshared"):
            nproc = int(value)
        elif name.lower() == "%mem":
            mem = value.strip()
    if nproc < 1:
        raise ValueError("%nproc/%nprocshared must be positive.")
    return dict(link0=link0, route=route, title=sections[1], charge=charge,
                multiplicity=multiplicity, atoms=atoms, geometry=np.array(coordinates),
                nproc=nproc, mem=mem, bottom="\n\n".join(sections[3:]))

##SOME CONSTANTS##############################################
EPSILON_0 = 8.854187817e-12  # F/m
HBAR_EV = 6.582119514e-16  # eV s
HBAR_J = 1.054571800e-34  # J s
MASS_E = 9.10938356e-31  # kg
LIGHT_SPEED = 299792458  # m/s
E_CHARGE = 1.60217662e-19  # C
BOLTZ_EV = 8.6173303e-5  # eV/K
AMU = 1.660539040e-27  # kg
###############################################################


##ERROR FUNCTION###############################################
def fatal_error(msg):
    print(msg)
    sys.exit()


###############################################################


##GETS FREQUENCIES AND REDUCED MASSES##########################
def pega_freq(freqlog):
    freqs, masses = [], []
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Frequencies -- " in line:
                line = line.split()
                for j in range(2, len(line)):
                    freqs.append(float(line[j]))
            elif "Red. masses --" in line:
                line = line.split()
                for j in range(3, len(line)):
                    masses.append(float(line[j]))
            elif "Thermochemistry" in line:
                break
    if len(freqs) == 0 or len(masses) == 0:
        fatal_error(
            "No frequencies or reduced masses found. Check your frequency file. Goodbye."
        )
    # conversion in angular frequency
    freqs = np.array(freqs) * (LIGHT_SPEED * 100 * 2 * np.pi)
    # conversion from amu to kg
    masses = np.asarray(masses) * AMU
    return freqs, masses


###############################################################


def pega_geom(freqlog):
    if os.path.splitext(str(freqlog))[1].lower() in (".com", ".gjf"):
        data = read_gaussian_input(freqlog)
        return data["geometry"], data["atoms"]
    geom = np.zeros((1, 3))
    atomos = []
    if ".log" in str(freqlog):
        busca = " orientation:"
        fetch = False
        with open(freqlog, "r",encoding='utf-8') as f:
            for line in f:
                if busca in line and "Dipole" not in line:
                    geom = np.zeros((1, 3))
                    atomos = []
                    fetch = True
                elif fetch:
                    line = line.split()
                    try:
                        float(line[0])
                        NG = []
                        for j in range(3, len(line)):
                            NG.append(float(line[j]))
                        atomos.append(line[1])
                        geom = np.vstack((geom, NG))
                    except (ValueError, IndexError):
                        if len(atomos) > 0:
                            fetch = False
    else:
        geom = np.zeros((1, 3))
        atomos = []
        with open(freqlog, "r",encoding='utf-8') as f:
            for line in f:
                line = line.split()
                try:
                    vetor = np.array([float(line[1]), float(line[2]), float(line[3])])
                    atomos.append(line[0])
                    geom = np.vstack((geom, vetor))
                except (IndexError, ValueError):
                    pass
    geom = geom[1:, :]
    atomos = [element_symbol(atom) for atom in atomos]
    return geom, atomos


###############################################################


##GETS NORMAL COORDINATES IN HIGH PRECISION####################
def pega_modosHP(G, freqlog):
    F, _ = pega_freq(freqlog)
    num_atoms = np.shape(G)[0]
    normal_modes = np.zeros((num_atoms, 3, len(F))) + np.nan
    mode = 0
    fetch = False
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Coord Atom Element:" in line:
                fetch = True
            elif fetch:
                line = line.split()
                if len(line) < 6 or "Harmonic" in line[0]:
                    fetch = False
                    mode += 5
                else:
                    coord = int(line[0]) - 1
                    atom = int(line[1]) - 1
                    for j in range(3, len(line)):
                        normal_modes[atom, coord, mode + j - 3] = float(line[j])
    return normal_modes


###############################################################


##GETS NORMAL COORDINATES IN REGULAR PRECISION#################
def pega_modosLP(G, freqlog):
    F, _ = pega_freq(freqlog)
    num_atoms = np.shape(G)[0]
    normal_modes = np.zeros((num_atoms, 3, len(F))) + np.nan
    mode = 0
    fetch = False
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Atom  AN      X      Y      Z" in line:
                fetch = True
            elif fetch:
                line = line.split()
                if len(line) < 5:
                    fetch = False
                    mode += 3
                else:
                    atom = int(line[0]) - 1
                    for j in range(2, len(line)):
                        normal_modes[atom, (j - 2) % 3, mode + (j - 2) // 3] = float(
                            line[j]
                        )
    return normal_modes


###############################################################


##DETECTS WHETHER HIGH PRECISION IS USED#######################
def pega_modos(geom, freqlog):
    x = "LP"
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Coord Atom Element:" in line:
                x = "HP"
                break
    if x == "LP":
        return pega_modosLP(geom, freqlog)
    else:
        return pega_modosHP(geom, freqlog)


###############################################################


##GETS INPUT PARAMS FROM LOG FILES#############################
def get_input_params(freqlog):
    nproc, mem, header = "", "", ""
    with open(freqlog, "r",encoding='utf-8') as f:
        search = False
        for line in f:
            if "%nproc" in line.lower():
                line = line.split("=")
                nproc = line[-1].replace("\n", "")
            elif "%mem" in line.lower():
                line = line.split("=")
                mem = line[-1].replace("\n", "")
            elif "#" in line and not search and header == "":
                search = True
                header += line.lstrip().replace("\n", "")
            elif search and "----------" not in line:
                header += line.lstrip().replace("\n", "")
            elif search and "----------" in line:
                search = False
                break
    return nproc, mem, header


###############################################################


##CHECKS THE FREQUENCY LOG'S LEVEL OF THEORY###################
def busca_input(freqlog):
    base = "lalala"
    exc = ""
    header = ""
    nproc = "4"
    mem = "1GB"
    scrf = ""
    nproc, mem, header = get_input_params(freqlog)

    if "TDA" in header.upper():
        exc = "tda"
        spec = "EMISPCT"
    elif "TD" in header.upper():
        exc = "td"
        spec = "EMISPCT"
    else:
        spec = "ABSSPCT"

    if "SCRF" in header.upper():
        new = header.split()
        for elem in new:
            if "SCRF" in elem:
                scrf = elem
                break

    header = header.split()
    base = ""
    for elem in header:
        if "/" in elem and "IOP" not in elem.upper():
            base += elem.replace("#", "")
        elif "IOP" in elem.upper() and ("108" in elem or "107" in elem):
            base += " " + elem
    return base, exc, nproc, mem, scrf, spec


###############################################################


##CHECKS FREQ FILES############################################
def double_check(freqlog):
    freqs, _ = pega_freq(freqlog)
    if np.any(freqs < 0):
        fatal_error(
            "Imaginary frequencies detected. Check your frequency file. Goodbye."
        )
    _, _, header = get_input_params(freqlog)
    header = header.lower()
    optfreqissue = False
    stationary = False
    if "opt" in header and "freq" in header and "iop(" in header:
        optfreqissue = True
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Stationary point found" in line:
                stationary = True
            elif all(s in line for s in ["Item", "Value", "Threshold", "Converged?"]):
                stationary = False
    if not stationary:
        print("*" * 50)
        print("WARNING: Non-optimized parameters detected in your frequency file.")
        print(
            "Even though the frequencies may be all real, your structure is still not fully optimized."
        )
        print("This may lead to inaccurate results.")
        if optfreqissue:
            print(
                "In your case, this may be due to running an opt freq calculation using a single input file and IOP options."
            )
            print(
                "Gaussian does not carry the IOP options to the frequency calculation when using a single input file."
            )
            print(
                "To avoid this issue, run the optimization and frequency calculations separately."
            )
        else:
            print("To learn more about this issue, check https://gaussian.com/faq3/ .")
        print("Proceed at your own risk.")
        print("*" * 50)
        print("\n")


###############################################################


##FETCHES REFRACTIVE INDEX#####################################
def get_nr():
    buscar = False
    coms = [
        file
        for file in os.listdir("Geometries")
        if "Geometr" in file and ".com" in file
    ]
    with open("Geometries/" + coms[0], "r",encoding='utf-8') as f:
        for line in f:
            if "SCRF" in line.upper():
                buscar = True
                break
    if buscar:
        logs = [
            file
            for file in os.listdir("Geometries")
            if "Geometr" in file and ".log" in file
        ]
        for log in logs:
            with open("Geometries/" + log, "r",encoding='utf-8') as f:
                for line in f:
                    if "Solvent" in line and "Eps" in line:
                        line = line.split()
                        nr = np.sqrt(float(line[6]))
                        return nr
    else:
        return 1


###############################################################


##FETCHES CHARGE AND MULTIPLICITY##############################
def get_cm(freqlog):
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "Charge" in line and "Multiplicity" in line:
                line = line.split()
                charge = line[2]
                mult = line[5]
                break
    return charge + " " + mult


###############################################################
