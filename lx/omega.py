#!/usr/bin/env python3
import os
import sys
import shutil
import numpy as np
import lx.tools
import lx.parser


##GENERATES NEUTRAL INPUT######################################
def gera_optcom(atomos, G, base, nproc, mem, omega, op):
    header = f"%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000) {op}\n\nTITLE\n\n0 1\n"
    header = (
        header.replace("JJJ", nproc)
        .replace("MEM", mem)
        .replace("BASE", base)
        .replace("MMMMM", omega)
    )
    file = "OPT_" + omega + "_.com"
    lx.tools.write_input(atomos, G, header, "", file)
    return file


###############################################################


##GENERATES ION INPUTS#########################################
def gera_ioncom(atomos, G, base, nproc, mem, omega):
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000)\n\nTITLE\n\n1 2\n"
    header = (
        header.replace("JJJ", nproc)
        .replace("MEM", mem)
        .replace("BASE", base)
        .replace("MMMMM", omega)
    )
    file1 = "pos_" + omega + "_.com"
    file2 = "neg_" + omega + "_.com"
    lx.tools.write_input(atomos, G, header, "", file1)
    header = "%nproc=JJJ\n%mem=MEM\n# BASE iop(3/108=MMMMM00000) iop(3/107=MMMMM00000) \n\nTITLE\n\n-1 2\n"
    header = (
        header.replace("JJJ", nproc)
        .replace("MEM", mem)
        .replace("BASE", base)
        .replace("MMMMM", omega)
    )
    lx.tools.write_input(atomos, G, header, "", file2)
    return [file1, file2]


###############################################################


##GETS ENERGY FROM LOG#########################################
def pega_energia(file):
    energias = []
    with open(file, "r",encoding='utf-8') as f:
        for line in f:
            if "SCF Done:" in line:
                line = line.split()
                energias.append(float(line[4]))
    return min(energias)


###############################################################


##MAP TO FLOAT WITH EXCEPTION##################################
def to_float(num):
    try:
        num = float(num)
    except ValueError:
        num = -np.inf
    return num


###############################################################


##GETS HOMO ENERGY FROM LOG####################################
def pega_homo(file):
    HOMOS = []
    with open(file, "r",encoding='utf-8') as f:
        for line in f:
            if "OPT" in file and "Optimized Parameters" in line:
                HOMOS = []
            if "occ. eigenvalues" in line:
                line = line.split()
                homos = line[4:]
                HOMOS.extend(homos)
    if len(HOMOS) == 0:
        with open(file, "r",encoding='utf-8') as f:
            for line in f:
                if "occ. eigenvalues" in line:
                    line = line.split()
                    homos = line[4:]
                    HOMOS.extend(homos)
    HOMOS = list(map(to_float, HOMOS))
    return np.max(HOMOS)


###############################################################


##RUNS CALCULATIONS############################################
def rodar_omega(atomos, geom, base, nproc, mem, omega, op, batch_file, gaussian, numjobs):
    omega = f"{omega:05.0f}"
    file = gera_optcom(atomos, geom, base, nproc, mem, omega, op)
    remover = []
    if op == "opt":
        the_watcher = lx.tools.Watcher('.',files=[file])
        the_watcher.run(batch_file, gaussian, 1)
        the_watcher.hold_watch()
        geom, atomos = lx.parser.pega_geom(file[:-3] + "log")
        files = gera_ioncom(atomos, geom, base, nproc, mem, omega)
        the_watcher = lx.tools.Watcher('.',files=files)
        the_watcher.run(batch_file, gaussian, min(numjobs,2))
        the_watcher.hold_watch()
    else:
        files = gera_ioncom(atomos, geom, base, nproc, mem, omega)
        the_watcher = lx.tools.Watcher('.',files=[file]+files)
        the_watcher.run(batch_file, gaussian, min(numjobs,3))
        the_watcher.hold_watch()

    logs = files + [file]
    for file in logs:
        remover.append(file)
        if "pos_" in file:
            cation = pega_energia(file[:-3] + "log")
        elif "neg_" in file:
            anion = pega_energia(file[:-3] + "log")
            homo_anion = pega_homo(file[:-3] + "log")
        elif "OPT_" in file:
            neutro = pega_energia(file[:-3] + "log")
            homo_neutro = pega_homo(file[:-3] + "log")

    try:
        os.mkdir("Logs")
    except FileExistsError:
        pass
    for file in remover:
        shutil.move(file, "Logs/" + file)
        shutil.move(file[:-3] + "log", "Logs/" + file[:-3] + "log")
    J = np.sqrt(
        ((homo_neutro + cation - neutro) ** 2 + (homo_anion + neutro - anion) ** 2)
    ) * (27.2114)
    return J


###############################################################


##WRITES LOG WITH RESULTS######################################
def write_tolog(omegas, Js, frase):
    with open("omega.lx", "w",encoding='utf-8') as f:
        f.write("#w(10^4 bohr^-1)    J(eV)\n")
        list1, list2 = zip(*sorted(zip(omegas, Js)))
        for i in range(len(list1)):
            f.write(f"{list1[i]:05.0f}               {list2[i]:.4f}\n")
        f.write(f"\n{frase} {list1[list2.index(min(list2))]:05.0f}\n")


###############################################################

def fetch_grad(x, f, i):
    """
    Estimate the gradient at x[i] using finite differences and Lagrange interpolation.

    Parameters:
    x (list or array): x values.
    f (list or array): Corresponding f(x) values.
    i (int): Index of the point to estimate the gradient.

    Returns:
    float: Estimated gradient at x[i].
    """
    if len(x) == 1:
        # Only one point; gradient is undefined, return 0
        return 0.0

    if len(x) == 2:
        # Two points; compute simple finite difference
        return (f[1] - f[0]) / (x[1] - x[0])

    if i == 0:
        # Forward finite difference for the first point
        return (f[1] - f[0]) / (x[1] - x[0])

    if i == len(x) - 1:
        # Backward finite difference for the last point
        return (f[-1] - f[-2]) / (x[-1] - x[-2])

    # General case: Use Lagrange interpolation for three points
    x0, x1, x2 = x[i - 1], x[i], x[i + 1]
    f0, f1, f2 = f[i - 1], f[i], f[i + 1]

    # Compute denominators for the Lagrange basis polynomials
    denom0 = (x0 - x1) * (x0 - x2)
    denom1 = (x1 - x0) * (x1 - x2)
    denom2 = (x2 - x0) * (x2 - x1)

    # First derivative (gradient) using corrected Lagrange formula
    l0p = (x1 - x2) / denom0
    l1p = (2 * x1 - x0 - x2) / denom1
    l2p = (x1 - x0) / denom2
    grad = f0 * l0p + f1 * l1p + f2 * l2p

    return grad

def main():
    geomlog = sys.argv[1]
    basis = sys.argv[2]
    nproc = sys.argv[3]
    mem = sys.argv[4]
    omega1 = sys.argv[5]
    passo = sys.argv[6]
    relax = sys.argv[7]
    script = sys.argv[8]
    gaussian = sys.argv[9]
    parallel = sys.argv[10]
    try:
        int(nproc)
        passo = float(passo) * 10000
        omega1 = float(omega1) * 10000
    except ValueError:
        lx.parser.fatal_error("nproc, omega and step must be numbers. Goodbye!")
    if relax.lower() == "y":
        op = "opt"
    elif relax.lower() == "n":
        op = ""
    else:
        lx.parser.fatal_error("It must be either y or n. Goodbye!")
    if parallel.lower() == "y":
        numjobs = 10
    else:
        numjobs = 1
    omegas, Js = [], []
    try:
        data = np.loadtxt("omega.lx", dtype=float)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        omegas = data[:, 0].tolist()
        Js = data[:, 1].tolist()
    except FileNotFoundError:
        pass    

    iteration = 0
    while iteration < 100:
        if omega1 in omegas:
            ind = omegas.index(omega1)
            J = Js[ind]
        else:
            try:
                #find existing omega closest to omega1
                d_omega = [abs(omega1 - om) for om in omegas]
                geom_log = omegas[d_omega.index(min(d_omega))]
                geom_log = f"Logs/OPT-{geom_log:3.0f}-.log"
                G, atomos = lx.parser.pega_geom(geom_log)
            except (ValueError, FileNotFoundError):
                G, atomos = lx.parser.pega_geom(geomlog)   
            J = rodar_omega(
                atomos, G, basis, nproc, mem, omega1, op, script, gaussian, numjobs)
            omegas.append(omega1)
            Js.append(J)
        omegas, Js = map(list, zip(*sorted(zip(omegas, Js))))
        #index of min J
        ind = Js.index(min(Js))
        omega1 = omegas[ind]
        grad = fetch_grad(omegas, Js, ind)
        
        if grad == 0:
            delta_omega = passo
        else:
            delta_omega = -1*Js[ind] / grad
        
        if ind == len(Js) - 1:
            max_omega = 5000
        else:
            max_omega = omegas[ind+1]
        if ind == 0:
            min_omega = 0
        else:
            min_omega = omegas[ind-1]
        
        omega1 += delta_omega
        omega1 = np.round(omega1, 0)
        if min_omega > omega1 or omega1 > max_omega or int(omega1) in omegas:
            left = (min_omega - omegas[ind])/2
            right = (max_omega - omegas[ind])/2
            omega1 =  omegas[ind] + max(left,right,key=abs)
        omega1 = int(omega1)
        if max_omega - min_omega <= 20:
            break


        write_tolog(omegas, Js, f"#Best value so far:")
        iteration += 1

    write_tolog(omegas, Js, "#Done! Optimized value:")
    menor = omegas[Js.index(min(Js))]
    log = f"Logs/OPT_{menor:05.0f}_.log"
    G, atomos = lx.parser.pega_geom(log)
    basis, _, nproc, mem, scrf, _ = lx.parser.busca_input(log)
    cm = lx.parser.get_cm(log)
    header = f"%nproc={nproc}\n%mem={mem}\n# {basis} {scrf}\n\nTITLE\n\n{cm}\n"
    lx.tools.write_input(atomos, G, header, "", "tuned_w.com")


if __name__ == "__main__":
    sys.exit(main())
