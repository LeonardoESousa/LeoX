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
    return J, geom, atomos


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


def main():
    geomlog = sys.argv[1]
    base = sys.argv[2]
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
    oms, jotas = [], []
    try:
        with open("omega.lx", "r",encoding='utf-8') as f:
            for line in f:
                line = line.split()
                if len(line) == 2 and "#" not in line:
                    om = float(line[0])
                    omegas.append(om)
                    Js.append(float(line[1]))
        menor = omegas[Js.index(min(Js))]
        G, atomos = lx.parser.pega_geom(f"Logs/OPT_{menor:05.0f}_.log")
    except FileNotFoundError:
        G, atomos = lx.parser.pega_geom(geomlog)

    while passo > 25:
        if omega1 in omegas:
            ind = omegas.index(omega1)
            J = Js[ind]
        else:
            J, G, atomos = rodar_omega(
                atomos, G, base, nproc, mem, omega1, op, script, gaussian, numjobs
            )
            omegas.append(omega1)
            Js.append(J)
        oms.append(omega1)
        jotas.append(J)
        try:
            if jotas[-1] - jotas[-2] > 0:
                passo = int(passo / 2)
                sign = -1 * np.sign(int(oms[-1]) - int(oms[-2]))
            else:
                sign = +1 * np.sign(int(oms[-1]) - int(oms[-2]))
            omega1 += sign * passo

        except IndexError:
            omega1 += passo
        write_tolog(omegas, Js, "#Best value so far:")

    write_tolog(omegas, Js, "#Done! Optimized value:")
    menor = omegas[Js.index(min(Js))]
    log = f"Logs/OPT_{menor:05.0f}_.log"
    G, atomos = lx.parser.pega_geom(log)
    base, _, nproc, mem, scrf, _ = lx.parser.busca_input(log)
    cm = lx.parser.get_cm(log)
    header = f"%nproc={nproc}\n%mem={mem}\n# {base} {scrf}\n\nTITLE\n\n{cm}\n"
    lx.tools.write_input(atomos, G, header, "", "tuned_w.com")


if __name__ == "__main__":
    sys.exit(main())
