#!/usr/bin/env python3
import os
import sys
import shutil
import numpy as np
import lx.tools
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


##SAMPLES GEOMETRIES###########################################
def make_geoms(freqlog, num_geoms, temp, header, bottom):
    lista = []
    counter = lx.tools.start_counter()
    _, atomos, structures = lx.tools.sample_geometries(
        freqlog, num_geoms, temp, 3000, warning=False
    )
    for n in range(np.shape(structures)[2]):
        final_geom = structures[:, :, n]
        lx.tools.write_input(
            atomos,
            final_geom,
            header.replace("UUUUU", str(n + 1)),
            bottom.replace("UUUUU", str(n + 1)),
            f"Geometry-{n+1+counter}-.com",
        )
        lista.append(f"Geometry-{n+1+counter}-.com")
    return lista


###############################################################


##GETS ENERGY FROM THE ORIGINAL FREQ LOG FILE##################
def get_energy_origin(freqlog):
    exc = 0
    with open(freqlog, "r",encoding='utf-8') as f:
        for line in f:
            if "SCF Done:" in line:
                line = line.split()
                scf = float(line[4]) * 27.2114
            elif "Total Energy," in line:
                line = line.split()
                exc = float(line[4]) * 27.2114
            elif "Rotational constants" in line:
                line = line.split()
                rot = [float(line[3]), float(line[4]), float(line[5])]
            elif "Normal termination" in line:
                if exc != 0:
                    scf = exc
                return scf, np.array([rot[0], rot[1], rot[2]])


###############################################################


##GETS ENERGIES FROM OPT LOG FILES#############################
def get_energies(folder, original_molecule):
    nums, scfs, rots = [], [], []
    files = [i for i in os.listdir(folder) if ".log" in i and "Geometry" in i]
    for file in files:
        exc = 0
        if np.array_equal(original_molecule, lx.tools.fingerprint(file, folder)):
            with open(folder + "/" + file, "r",encoding='utf-8') as f:
                num = float(file.split("-")[1])
                for line in f:
                    if "SCF Done:" in line:
                        line = line.split()
                        scf = float(line[4]) * 27.2114
                    elif "Total Energy," in line:
                        line = line.split()
                        exc = float(line[4]) * 27.2114
                    elif "Rotational constants" in line:
                        line = line.split()
                        rot = [float(line[3]), float(line[4]), float(line[5])]
                    elif "Normal termination" in line:
                        if exc != 0:
                            scf = exc
                        scfs.append(scf)
                        nums.append(num)
                        rots.append(rot)

    for file in files:
        try:
            shutil.move(file, "Geometries/" + file)
            shutil.move(file[:-3] + "com", "Geometries/" + file[:-3] + "com")
        except FileNotFoundError:
            pass
    return np.array(nums), np.array(scfs), np.array(rots)


###############################################################


def measure(vec1, vec2, cr):
    vec1 = np.array(vec1)
    vec2 = np.array(vec2)
    cr = np.array(cr)
    dist = max(abs(vec1 - vec2) - cr)
    distance = np.heaviside(dist, 0)
    return distance


class Conformation:
    def __init__(self, rot, energy, identity, num) -> None:
        self.rot = rot[np.newaxis, :]
        self.energy = [energy]
        self.identity = identity
        self.std = rot / 1000
        self.num = [num]

    def add_rot(self, rot, num, energy):
        self.rot = np.vstack((self.rot, rot[np.newaxis, :]))
        newstd = np.std(self.rot, axis=0)
        # get higher std
        self.std = np.where(newstd > self.std, newstd, self.std)
        self.num.append(num)
        self.energy.append(energy)

    def merge_rot(self, rot):
        self.rot = np.vstack((self.rot, rot))
        newstd = np.std(self.rot, axis=0)
        # get higher std
        self.std = np.where(newstd > self.std, newstd, self.std)

    def get_avg(self):
        return np.mean(self.rot, axis=0)


def internal_comparison(conformations):
    num_conformations = len(conformations)
    remove = set()  # Use a set to avoid duplicates

    for i in range(num_conformations - 1):
        if i in remove:
            continue  # Skip already marked for removal
        for j in range(i + 1, num_conformations):
            if j in remove:
                continue  # Skip already marked for removal
            distance = measure(
                conformations[i].get_avg(),
                conformations[j].get_avg(),
                conformations[j].std + conformations[i].std,
            )
            if distance == 0:
                conformations[i].num += conformations[j].num
                conformations[i].energy += conformations[j].energy
                conformations[i].merge_rot(conformations[j].rot)
                remove.add(j)
                
    # Remove elements from list whose index is in remove
    conformations = [c for idx, c in enumerate(conformations) if idx not in remove]
    return conformations


def classify(conformations, folder):
    nums, scfs, rots = get_energies(folder, conformations[0].identity)
    new_conformations = []

    total_geometries = rots.shape[0]
    matched_count = 0
    new_conformation_count = 0

    for i in range(total_geometries):
        matched = False
        for conformation in conformations:
            distance = measure(rots[i, :], conformation.get_avg(), conformation.std)
            if distance == 0:
                conformation.add_rot(rots[i, :], nums[i], scfs[i])
                matched = True
                matched_count += 1
                break

        if not matched:
            new_conformations.append(Conformation(rots[i, :], scfs[i], conformations[0].identity, nums[i]))
            new_conformation_count += 1

    conformations.extend(new_conformations)

    if len(conformations) > 1:
        conformations = internal_comparison(conformations)
    return conformations



def write_report(conformations, rounding, total_rounds, temp):
    engs = []
    for conformation in conformations:
        engs.append(np.mean(conformation.energy))
    engs = np.array(engs)
    argsort = np.argsort(engs)
    engs = engs[argsort]
    # sort conformations as list
    conformations = [conformations[i] for i in argsort]
    deltae = engs - min(engs)
    probs = np.exp(-(deltae) / (0.026))
    probs = probs / sum(probs)

    with open("conformation.lx", "w",encoding='utf-8') as f:
        f.write(
            f"{'#Group':<6}  {'Energy(eV)':<10}  {'DeltaE(eV)':<10}  {'Prob@300K(%)':<12}  {'Rot1':<10}  {'Rot2':<10}  {'Rot3':<10}  {'Std1':<10}  {'Std2':<10}  {'Std3':<10}  {'Number':<6}  {'Last':<6}\n"
        )
        for i, conformation in enumerate(conformations):
            rot = conformation.get_avg()
            std = conformation.std
            last = conformation.num[-1]
            total = len(conformation.num)
            f.write(
                f"{i+1:<6}  {engs[i]:<10.3f}  {deltae[i]:<10.3f}  {100*probs[i]:<12.1f}  {rot[0]:<10.7f}  {rot[1]:<10.7f}  {rot[2]:<10.7f}  {std[0]:<10.7f}  {std[1]:<10.7f}  {std[2]:<10.7f}  {total:<6.0f}  {last:<6.0f}\n"
            )
        f.write(f"\n#Round {rounding}/{total_rounds} Temperature: {temp} K")


##RUNS FREQ CALCULATION FOR NEW CONFORMATION###################
def rodar_freq(origin, nproc, mem, base, cm, batch_file, gaussian):
    geomlog = f"Geometries/Geometry-{origin:.0f}-.log"
    geom, atomos = lx.parser.pega_geom(geomlog)
    header = f"%nproc={nproc}\n%mem={mem}\n# freq=(noraman) nosymm  {base} \n\nTITLE\n\n{cm}\n"
    file = f"Freq-{origin:.0f}-.com"
    lx.tools.write_input(atomos, geom, header, "", file)
    the_watcher = lx.tools.Watcher('.',files=[file])
    the_watcher.run(batch_file, gaussian, 1)
    the_watcher.hold_watch()
    #lx.tools.rodar_lista([file], batch_file, gaussian, "conformation.lx", 1)
    log = file[:-3] + "log"
    with open(log, "r",encoding='utf-8') as f:
        for line in f:
            if "Normal termination" in line:
                return log
            elif "Error termination" in line:
                return None


###############################################################


def classify_only():
    files = [i for i in os.listdir(".") if "Geometry-" in i and ".log" in i]
    conformations = []
    for file in files:
        num = int(file.split("-")[1])
        scf, rot = get_energy_origin(file)
        conformations.append(
            Conformation(rot, scf, lx.tools.fingerprint(file, "."), num)
        )
    conformations = internal_comparison(conformations)
    write_report(conformations, 0, 0, 0)


def main():
    freqlog = sys.argv[1]
    base = sys.argv[2]
    nproc = sys.argv[3]
    mem = sys.argv[4]
    temperature = float(sys.argv[5])
    delta_t = float(sys.argv[6])
    num_geoms = int(sys.argv[7])
    rounds = int(sys.argv[8])
    numjobs = int(sys.argv[9])
    script = sys.argv[10]
    gaussian = sys.argv[11]
    temp_0 = temperature
    freq0 = freqlog
    if "td" in base.lower():
        opt = "=loose"
    else:
        opt = ""

    try:
        os.mkdir("Geometries")
    except FileExistsError:
        pass
    cm = lx.parser.get_cm(freqlog)
    header = (
        f"%nproc={nproc}\n%mem={mem}\n# opt{opt} nosymm  {base} \n\nTitle\n\n{cm}\n"
    )
    scf, rot = get_energy_origin(freqlog)
    conformations = [Conformation(rot, scf, lx.tools.fingerprint(freqlog, "."), 0)]
    files = [i for i in os.listdir("Geometries") if "Geometry" in i and ".log" in i]
    if len(files) > 0:
        conformations = classify(conformations, "Geometries")
        write_report(conformations, 0, rounds, temp_0)
    else:
        pass

    groups = len(conformations)
    for i in range(rounds):
        lista = make_geoms(freqlog, num_geoms, temp_0, header, "")
        the_watcher = lx.tools.Watcher('.',files=lista)
        the_watcher.run(script, gaussian, numjobs)
        the_watcher.hold_watch()
        #lx.tools.rodar_lista(lista, script, gaussian, "conformation.lx", numjobs)
        conformations = classify(conformations, ".")
        write_report(conformations, i + 1, rounds, temp_0)

        if len(conformations) != groups:
            log = rodar_freq(
                conformations[-1].num[-1], nproc, mem, base, cm, script, gaussian
            )
            if log is not None:
                freqlog = log
                temp_0 = temperature
            groups = len(conformations)
        else:
            temp_0 += delta_t

    with open("conformation.lx", "a",encoding='utf-8') as f:
        f.write("\n#Search concluded!")

    try:
        os.mkdir("Conformers")
    except FileExistsError:
        pass

    for i, conformation in enumerate(conformations):
        numero = conformation.num[-1]
        if numero == 0:
            freqlog = freq0
        else:
            freqlog = f"Geometries/Geometry-{numero:.0f}-.log"
        _, _, nproc, mem, scrf, _ = lx.parser.busca_input(freqlog)
        cm = lx.parser.get_cm(freqlog)
        header = f"%nproc={nproc}\n%mem={mem}\n%chk=Group_{i+1}_.chk\n# pm6 {scrf} opt nosymm\n\nTITLE\n\n{cm}\n"
        geom, atomos = lx.parser.pega_geom(freqlog)
        lx.tools.write_input(
            atomos, geom, header, "", f"Conformers/Geometry-{i+1}-.com"
        )


if __name__ == "__main__":
    sys.exit(main())
