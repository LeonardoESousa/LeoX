#!/usr/bin/env python3
"""xTB/CREST sampling followed by Gaussian refinement and CREGEN sorting."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from tempfile import TemporaryDirectory

import numpy as np
import lx.parser
import lx.tools

HARTREE_EV = 27.211386245988
TEMPERATURE = 300.0
PRUNING_VERSION = 4


def read_xyz(filename):
    """Read an XYZ ensemble, retaining the energy comment in each block."""
    lines = Path(filename).read_text(encoding="utf-8").splitlines()
    structures, index = [], 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        count = int(lines[index])
        if count < 1 or index + count + 2 > len(lines):
            raise ValueError("Incomplete XYZ block in " + str(filename))
        atoms, geometry = [], []
        for line in lines[index + 2:index + count + 2]:
            fields = line.split()
            if len(fields) != 4:
                raise ValueError("Invalid XYZ coordinate in " + str(filename))
            atoms.append(lx.parser.element_symbol(fields[0]))
            geometry.append([float(value) for value in fields[1:]])
        if not np.isfinite(geometry).all():
            raise ValueError("Nonfinite XYZ coordinates.")
        structures.append(dict(atoms=atoms, geometry=np.array(geometry), comment=lines[index + 1]))
        index += count + 2
    if not structures:
        raise ValueError("Empty ensemble: " + str(filename))
    return structures


def write_xyz(filename, structures):
    with open(filename, "w", encoding="utf-8") as handle:
        for structure in structures:
            handle.write(f"{len(structure['atoms'])}\n{structure.get('comment', '')}\n")
            for atom, xyz in zip(structure["atoms"], structure["geometry"]):
                handle.write(f"{atom} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f}\n")


def run_command(command, folder, logfile):
    with open(Path(folder) / logfile, "w", encoding="utf-8") as output:
        subprocess.run(command, cwd=folder, stdout=output, stderr=subprocess.STDOUT, check=True)


def sbatch_command(template, script, folder, nproc, mem):
    match = re.fullmatch(r"(\d+)\s*([KMGT]?)([BW]?)", mem, re.I)
    if not match:
        raise ValueError("Unsupported Gaussian memory specification: " + mem)
    amount, unit, kind = match.groups()
    amount = int(amount)
    unit, kind = unit.upper(), kind.upper()
    if not unit:
        # Gaussian's unitless memory is in 8-byte words; SLURM defaults to MB.
        slurm_mem = str(max(1, (amount * 8 + 1048575) // 1048576)) + "M"
    else:
        slurm_mem = str(amount * (8 if kind == "W" else 1)) + unit
    return ["sbatch", "--wait", "--parsable", "--nodes=1", "--ntasks=1",
            "--ntasks-per-node=1", f"--cpus-per-task={nproc}", f"--mem={slurm_mem}",
            f"--chdir={Path(folder).resolve()}", str(Path(template).resolve()), str(Path(script).resolve())]


def run_slurm(template, script, folder, nproc, mem):
    """Wait for completion and propagate submission or job failure."""
    command = sbatch_command(template, script, folder, nproc, mem)
    with open(Path(folder) / "slurm_submission.log", "w", encoding="utf-8") as output:
        process = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT)
        while process.poll() is None:
            time.sleep(20)
        if process.returncode:
            raise RuntimeError(f"SLURM job failed ({process.returncode}); see {folder}/slurm_submission.log")


def route_tokens(route):
    """Keep parenthesized Gaussian option lists together, including spaces."""
    tokens, token, depth = [], "", 0
    route = re.sub(r"^\s*#[pPnNtT]?\s*", "", route)
    route = re.sub(r"\s*=\s*", "=", route)
    for char in route:
        if char.isspace() and depth == 0:
            if token:
                tokens.append(token)
                token = ""
        else:
            token += char
            depth += (char == "(") - (char == ")")
    if token:
        tokens.append(token)
    if depth:
        raise ValueError("Unbalanced Gaussian route options.")
    return tokens


def gaussian_routes(template):
    common, opt, freq_option, connected = [], "opt", "freq=noraman", False
    for token in route_tokens(template["route"]):
        key = re.split(r"[=(]", token.lower(), maxsplit=1)[0]
        if key == "opt":
            if re.search(r"\b(ts|qst2|qst3|restart|modredundant|addgic|readfreeze)\b", token, re.I):
                raise ValueError("Use an unconstrained minimum optimization template, not TS/QST/restart/read constraints.")
            opt = token
        elif key == "freq":
            if re.search(r"\b(readfc|readisotopes|fc|fcht|anharmonic)\b", token, re.I):
                raise ValueError("Use a fresh harmonic frequency calculation, not a checkpoint/property frequency job.")
            freq_option = token
        elif key == "temperature":
            continue
        elif key == "geom":
            if token.lower() != "geom=connectivity":
                raise ValueError("Only explicit Cartesian geometry (optionally geom=connectivity) is supported.")
            connected = True
        elif key == "units" and token.lower() not in ("units=angstroms", "units=angstrom"):
            raise ValueError("Supply Cartesian coordinates in Angstroms.")
        else:
            common.append(token)
    route = " ".join(common)
    optimization = f"#p {route} {opt}" + (" geom=connectivity" if connected else "")
    frequency_route = " ".join(token for token in common if re.split(r"[=(]", token.lower(), maxsplit=1)[0] != "guess")
    frequency = f"#p {frequency_route} {freq_option} temperature=300 geom=allcheck guess=read"
    tail = template["bottom"]
    if connected:
        # AllCheck does not read the explicit connectivity block again.
        blocks = re.split(r"\n\s*\n", tail, maxsplit=1)
        tail = blocks[1] if len(blocks) == 2 else ""
    return optimization, frequency, tail


def make_gaussian_inputs(template, structures, folder):
    optimization, frequency, frequency_tail = gaussian_routes(template)
    files = []
    for index, structure in enumerate(structures, 1):
        if structure["atoms"] != template["atoms"]:
            raise ValueError("CREST changed the atom sequence; cannot reuse the Gaussian template.")
        index = structure.get("crest_index", index)
        name = f"Geometry-{index}-.com"
        if (Path(folder) / name).is_file():
            files.append(name)
            continue
        checkpoint = f"conformer_{index}.chk"
        # Use a distinct checkpoint/scratch file for every concurrent calculation.
        link0 = [line for line in template["link0"]
                 if line.partition("=")[0].lower() not in ("%chk", "%rwf", "%int", "%d2e", "%nproc", "%nprocshared", "%mem")]
        link0 += [f"%nprocshared={template['nproc']}", f"%mem={template['mem']}", f"%chk={checkpoint}"]
        cm = f"{template['charge']} {template['multiplicity']}"
        header = "\n".join(link0) + f"\n{optimization}\n\n{template['title']} - conformer {index}\n\n{cm}\n"
        bottom = template["bottom"] + "\n\n--Link1--\n"
        bottom += f"%nprocshared={template['nproc']}\n%mem={template['mem']}\n%chk={checkpoint}\n"
        bottom += frequency + "\n\n" + frequency_tail + "\n\n"
        lx.tools.write_input(structure["atoms"], structure["geometry"], header, bottom, str(Path(folder) / name))
        files.append(name)
    return files


class ConformerWatcher(lx.tools.Watcher):
    """Use Watcher for Gaussian completion, also detecting killed SLURM jobs."""
    def __init__(self, files):
        super().__init__(".", files=files, counter=2)
        self.status_files = {name[:-4]: f"cmd_{i}_.sh.status" for i, name in enumerate(files)}

    def check(self):
        super().check()
        for name in self.files.copy():
            # sbatch --wait finished but two Gaussian terminations never arrived.
            if Path(self.status_files[name]).exists():
                self.error.append(name)
                self.files.remove(name)


def run_gaussian_jobs(template, batch, files, folder, max_jobs):
    # Watcher calls `bash wrapper script`. A background waiting sbatch process
    # preserves concurrency while recording timeouts/OOM/cancelled jobs.
    folder = Path(folder).resolve()
    command = sbatch_command(batch, "unused", folder, template["nproc"], template["mem"])[:-1]
    submission = " ".join(shlex.quote(arg) for arg in command) + ' "$1"'
    wrapper = folder / "submit_gaussian.sh"
    wrapper.write_text('#!/bin/bash\n(\n' + submission + ' > "$1.slurm.log" 2>&1\n'
                       'status=$?\nprintf "%s\\n" "$status" > "$1.status"\n) < /dev/null &\n', encoding="utf-8")
    (folder / "limit.lx").write_text(str(max_jobs), encoding="utf-8")
    # Old logs and status files must not make Watcher finish a resubmitted job early.
    for index, name in enumerate(files):
        (folder / f"cmd_{index}_.sh.status").unlink(missing_ok=True)
        (folder / Path(name).with_suffix(".log")).unlink(missing_ok=True)
    previous = Path.cwd()
    os.chdir(folder)
    try:
        watcher = ConformerWatcher(files)
        watcher.run(str(wrapper), template["gaussian"], 1)
        watcher.hold_watch()
        return watcher.error
    finally:
        os.chdir(previous)


def frequency_is_stationary(filename):
    """Check the last opt/freq convergence check, not just real frequencies."""
    text = Path(filename).read_text(encoding="utf-8", errors="replace")
    if "Error termination" in text or text.count("Normal termination") != 2:
        return False
    freq = text.split("Normal termination")[1]
    values = [float(value) for line in freq.splitlines() if "Frequencies --" in line
              for value in line.split("--", 1)[1].split()]
    if not values or not np.isfinite(values).all() or any(value < 0 for value in values):
        return False
    stationary = False
    for line in text.splitlines():
        if "Stationary point found" in line or "Optimization completed" in line:
            stationary = True
        elif all(word in line for word in ("Item", "Value", "Threshold", "Converged?")):
            stationary = False
    return stationary


def save_json(filename, data):
    """Commit stage metadata atomically so an interruption cannot truncate it."""
    filename = Path(filename)
    temporary = filename.with_name(filename.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2), encoding="utf-8")
    temporary.replace(filename)


def file_digest(filename):
    import hashlib
    digest = hashlib.sha256()
    with open(filename, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gaussian_finished(filename):
    """Completion is distinct from passing the stationary-point/minimum check."""
    filename = Path(filename)
    if not filename.is_file():
        return False
    text = filename.read_text(encoding="utf-8", errors="replace")
    return ("Error termination" not in text and text.count("Normal termination") == 2
            and "Frequencies --" in text.split("Normal termination")[1]
            and "Sum of electronic and thermal Free Energies" in text.split("Normal termination")[1])


def sampling_finished(folder, output, logfile, termination, atoms, single=False):
    try:
        structures = read_xyz(folder / output)
        text = (folder / logfile).read_text(encoding="utf-8", errors="replace").lower()
        return (termination.lower() in text
                and (not single or (len(structures) == 1 and "geometry optimization converged" in text))
                and all(item["atoms"] == atoms for item in structures))
    except (OSError, ValueError, IndexError):
        return False


def remove_output(path):
    """Remove a generated output without following directory symlinks."""
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def restart_stage(folder):
    """Discard incomplete stage outputs and restart in the same folder."""
    for path in folder.iterdir():
        remove_output(path)


def classification_signature(folder, files, crest, rthr, ethr, merge_mirrors=True):
    files = sorted(set(files) | {path.with_suffix(".com").name
                                for path in (folder / "Geometries").glob("Geometry-*.log")})
    return dict(logs={name: file_digest(folder / "Geometries" / Path(name).with_suffix(".log"))
                      if (folder / "Geometries" / Path(name).with_suffix(".log")).is_file()
                      else None for name in files},
                crest=str(crest), rthr=rthr, ethr=ethr, temperature=TEMPERATURE,
                pruning=PRUNING_VERSION, merge_mirrors=merge_mirrors,
                sampling_groups=file_digest(folder / "CREST" / "symmetry_groups.json")
                if (folder / "CREST" / "symmetry_groups.json").exists() else None)


def cached_classification(folder, signature):
    marker = folder / "CREGEN" / "classification.json"
    try:
        state = json.loads(marker.read_text(encoding="utf-8"))
        if state["signature"] != signature:
            return None
        for name, digest in state["outputs"].items():
            if file_digest(folder / name) != digest:
                return None
        if set(state["outputs"]) != {"conformation.csv", "conformers_manifest.csv", "conformers_unique.xyz"}:
            return None
        read_xyz(folder / "conformers_unique.xyz")
        with open(folder / "conformation.csv", newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        with open(folder / "conformers_manifest.csv", newline="", encoding="utf-8") as handle:
            metadata = {row["Gaussian_log"]: row for row in csv.DictReader(handle)}
        results = [gaussian_result(folder / "Geometries" / row["Gaussian_log"]) for row in rows]
        for item, row in zip(results, rows):
            item["multiplicity"] = int(row["Multiplicity"])
            details = metadata[row["Gaussian_log"]]
            item["enantiomers"] = details["Enantiomers"] == "yes"
            item["partner_log"] = details["Partner_Gaussian_log"]
            item["minima"] = [dict(item)]
            if item["multiplicity"] == 2:
                partner = item["partner_log"]
                item["minima"].append(gaussian_result(folder / "Geometries" / partner)
                                      if partner != "inferred mirror partner" else item["minima"][0])
        return results or None
    except (OSError, ValueError, KeyError, TypeError):
        return None


def retry_frequency_checks(template, batch, files, folder, max_jobs):
    """One completed retry per conformer; resume interrupted retries in place."""
    from shutil import copy2

    folder = Path(folder).resolve()
    statefile = folder / "frequency_retries.json"
    state = json.loads(statefile.read_text(encoding="utf-8")) if statefile.exists() else {}
    originals = folder / "Retry" / "Originals"
    _, frequency, tail = gaussian_routes(template)
    common = [token for token in route_tokens(frequency)
              if re.split(r"[=(]", token.lower(), maxsplit=1)[0]
              not in ("freq", "temperature", "geom", "guess")]
    optimization = "#p " + " ".join(common) + " opt=readfc guess=read geom=allcheck"
    retries = []
    for name in files:
        original = folder / Path(name).with_suffix(".log")
        checkpoint = f"conformer_{int(Path(name).stem.split('-')[1])}.chk"
        entry = state.get(name)
        # Recognize second attempts made by the previous incremental patch.
        legacy = folder / "Retry" / name
        if entry is None and legacy.is_file() and gaussian_finished(legacy.with_suffix(".log")):
            originals.mkdir(parents=True, exist_ok=True)
            try:
                gaussian_result(legacy.with_suffix(".log"))
                for path in (original, folder / checkpoint):
                    if path.is_file() and not (originals / path.name).exists():
                        copy2(path, originals / path.name)
                copy2(legacy.with_suffix(".log"), original)
                if legacy.with_name(checkpoint).is_file():
                    copy2(legacy.with_name(checkpoint), folder / checkpoint)
            except (OSError, ValueError):
                pass  # Preserve the first attempt if the old retry was unusable.
            (folder / name).write_text(re.sub(r"(?im)^%oldchk=.*\n", "", legacy.read_text()), encoding="utf-8")
            state[name] = dict(phase="done")
            save_json(statefile, state)
            entry = state[name]
        if entry and entry["phase"] == "done":
            if gaussian_finished(original) and not frequency_is_stationary(original):
                print(f"WARNING: {original.name} still fails the frequency check; "
                      "its second attempt has already finished. Continuing.", flush=True)
            continue
        if entry is None:
            if not gaussian_finished(original) or frequency_is_stationary(original):
                continue
            if not (folder / checkpoint).is_file():
                print(f"WARNING: {original.name} failed the frequency check; "
                      "checkpoint missing, continuing without a retry.", flush=True)
                continue
            originals.mkdir(parents=True, exist_ok=True)
            for path in (original, folder / checkpoint):
                target = originals / path.name
                if not target.exists():
                    copy2(path, target)
            link0 = (f"%nprocshared={template['nproc']}\n%mem={template['mem']}\n"
                     f"%chk={checkpoint}\n")
            retry_input = link0 + optimization + "\n\n" + tail + "\n\n--Link1--\n"
            retry_input += link0 + frequency + "\n\n" + tail + "\n\n"
            # Save the intended input before modifying any job files.
            state[name] = dict(phase="prepared", original_log=file_digest(original), input=retry_input)
            save_json(statefile, state)
            entry = state[name]
        target = folder / name
        temporary = target.with_name(target.name + ".tmp")
        temporary.write_text(entry["input"], encoding="utf-8")
        temporary.replace(target)
        # A crash between preparing the retry and removing the old log is harmless.
        if original.is_file() and file_digest(original) == entry["original_log"]:
            original.unlink()
        if not gaussian_finished(original):
            # Restore the frequency Hessian if an interrupted retry changed the checkpoint.
            copy2(originals / checkpoint, folder / checkpoint)
            retries.append(name)
    if retries:
        print(f"Submitting {len(retries)} Gaussian frequency retries with Opt=ReadFC...", flush=True)
        run_gaussian_jobs(template, batch, retries, folder, max_jobs)
    for name, entry in state.items():
        if entry["phase"] != "prepared":
            continue
        original = folder / Path(name).with_suffix(".log")
        if not gaussian_finished(original):
            text = original.read_text(encoding="utf-8", errors="replace") if original.exists() else ""
            if "Error termination" not in text:
                print(f"WARNING: Retry for {name} unfinished; it will resume on the next launch.", flush=True)
                continue
            print(f"WARNING: Retry for {name} ended with an error; keeping the first attempt.", flush=True)
            checkpoint = f"conformer_{int(Path(name).stem.split('-')[1])}.chk"
            copy2(originals / original.name, original)
            copy2(originals / checkpoint, folder / checkpoint)
            entry["phase"] = "done"
            save_json(statefile, state)
            continue
        try:
            gaussian_result(original)
        except (OSError, ValueError) as error:
            print(f"WARNING: Retry for {name} unusable ({error}); keeping the first attempt.", flush=True)
            checkpoint = f"conformer_{int(Path(name).stem.split('-')[1])}.chk"
            copy2(originals / original.name, original)
            copy2(originals / checkpoint, folder / checkpoint)
        if not frequency_is_stationary(original):
            print(f"WARNING: {original.name} still fails the frequency convergence "
                  "check after one retry; continuing.", flush=True)
        entry["phase"] = "done"
        save_json(statefile, state)
    # These are recovery files, needed only while a second attempt is unfinished.
    if all(entry["phase"] == "done" for entry in state.values()):
        remove_output(folder / "Retry")


def gaussian_result(filename):
    """Read completed opt/freq results, in Hartree and Kelvin."""
    text = Path(filename).read_text(encoding="utf-8", errors="replace")
    if "Error termination" in text or text.count("Normal termination") != 2:
        raise ValueError("Gaussian opt/freq did not finish normally twice.")
    jobs = text.split("Normal termination")
    opt, freq = jobs[0], jobs[1]
    if "Optimization completed" not in opt and "Stationary point found" not in opt:
        raise ValueError("Gaussian optimization did not reach a stationary point.")
    frequencies = []
    for line in freq.splitlines():
        if "Frequencies --" in line:
            frequencies.extend(float(value) for value in line.split("--", 1)[1].split())
    if not frequencies or any(value < 0 for value in frequencies):
        raise ValueError("Missing frequencies or imaginary frequencies; not a verified minimum.")
    energies = re.findall(r"SCF Done:.*?=\s*([-+\d.DEde]+)", opt)
    excited = re.findall(r"Total Energy,\s*E\([^)]*\)\s*=\s*([-+\d.DEde]+)", opt)
    gibbs = re.findall(r"Sum of electronic and thermal Free Energies\s*=\s*([-+\d.DEde]+)", freq)
    corrections = re.findall(r"Thermal correction to Gibbs Free Energy\s*=\s*([-+\d.DEde]+)", freq)
    temperatures = re.findall(r"Temperature\s+([\d.]+)\s+Kelvin", freq)
    if not gibbs or not corrections or not temperatures:
        raise ValueError("Missing electronic energy, Gibbs energy, or thermochemistry temperature.")
    if abs(float(temperatures[-1]) - TEMPERATURE) > 0.01:
        raise ValueError("Thermochemistry was not evaluated at 300 K.")
    convert = lambda value: float(value.replace("D", "E").replace("d", "e"))
    free_energy = convert(gibbs[-1])
    # Thermochemistry supplies the correct reference energy even for correlated
    # methods. Last displaced SCF energies in numerical Hessians are unsuitable.
    electronic = free_energy - convert(corrections[-1])
    candidate = excited[-1] if excited else (energies[-1] if energies else None)
    if candidate is not None and abs(convert(candidate) - electronic) < 2e-6:
        electronic = convert(candidate)
    if not np.isfinite([electronic, free_energy]).all():
        raise ValueError("Nonfinite Gaussian energy.")
    # Frequency calculations may end with displaced geometries. Use the last
    # orientation of the optimization job rather than the last in the whole log.
    with TemporaryDirectory() as temporary:
        opt_log = Path(temporary) / "opt.log"
        opt_log.write_text(opt, encoding="utf-8")
        geometry, atoms = lx.parser.pega_geom(str(opt_log))
    if not len(atoms) or not np.isfinite(geometry).all():
        raise ValueError("Missing final Gaussian geometry.")
    return dict(atoms=atoms, geometry=geometry, energy=electronic, gibbs=free_energy,
                comment=f"{electronic:.12f}", source=str(Path(filename).resolve()))


def aligned_rmsd(geometry, reference):
    """Compare geometries after translation and proper rotation (no reflection).

    This identifies the source of a CREGEN output, not a new clustering rule.
    CREGEN's own thresholds still determine which conformers are retained.
    """
    geometry = np.asarray(geometry, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if geometry.shape != reference.shape or geometry.size == 0:
        return float("inf")
    centered = geometry - geometry.mean(axis=0)
    target = reference - reference.mean(axis=0)
    left, _, right = np.linalg.svd(centered.T @ target)
    if np.linalg.det(left @ right) < 0:
        left[:, -1] *= -1
    difference = centered @ (left @ right) - target
    return float(np.sqrt(np.mean(np.sum(difference ** 2, axis=1))))


def symmetry_rmsd(first, second, inversion=False):
    """iRMSD with an explicit inversion policy and connectivity verification."""
    import irmsd
    if sorted(first["atoms"]) != sorted(second["atoms"]):
        return float("inf")
    # The 0.1.2 backend can crash on one/two-atom principal-axis degeneracies.
    # Their rotation-invariant distance is analytic; neither can be chiral.
    if len(first["atoms"]) <= 2:
        if len(first["atoms"]) == 1:
            return 0.0
        distances = [np.linalg.norm(item["geometry"][1] - item["geometry"][0])
                     for item in (first, second)]
        return float(abs(distances[0] - distances[1]) / 2)
    molecules = [irmsd.Molecule(item["atoms"], item["geometry"]) for item in (first, second)]
    value, left, right = irmsd.get_irmsd_molecule(*molecules, iinversion=1 if inversion else 2)
    if not np.isfinite(value):
        raise ValueError("iRMSD returned a nonfinite distance.")
    # Canonical ranks alone need not guarantee a graph-preserving assignment.
    # Verify the connectivity in the aligned, reordered coordinates as well.
    if not np.array_equal(lx.tools.adjacency(left.positions, left.symbols),
                          lx.tools.adjacency(right.positions, right.symbols)):
        return float("inf")
    return float(value)


def prune_symmetry(structures, rthr, ethr, merge_mirrors=True):
    """Group duplicates and observed mirror partners; never count sampling repeats."""
    def energy(item):
        if "energy" in item:
            return item["energy"]
        try:
            return float(item.get("comment", "").split()[0])
        except (ValueError, IndexError):
            return 0.0
    groups = []
    for item in sorted(structures, key=energy):
        for group in groups:
            representative = group["representative"]
            proper = symmetry_rmsd(representative, item)
            mirrored = proper > rthr and merge_mirrors and symmetry_rmsd(representative, item, True) <= rthr
            if proper > rthr and not mirrored:
                continue
            if abs(energy(item) - energy(representative)) * 627.509474 > ethr + 1e-10:
                print("WARNING: Similar geometries have inconsistent electronic energies; keeping both.", flush=True)
                continue
            group["members"].append((item, "enantiomer" if mirrored else "duplicate"))
            # One minimum per handedness, irrespective of how often it was sampled.
            if mirrored and len(group["minima"]) == 1:
                group["minima"].append(item)
            break
        else:
            groups.append(dict(representative=item, members=[(item, "representative")], minima=[item]))
    return groups


def sampling_prune(structures, folder, rthr, ethr, merge_mirrors):
    for index, item in enumerate(structures, 1):
        item["crest_index"] = index
    groups = prune_symmetry(structures, rthr, ethr, merge_mirrors)
    metadata = []
    for group in groups:
        representative = group["representative"]
        metadata.append(dict(gaussian_log=f"Geometry-{representative['crest_index']}-.log",
                             enantiomers=len(group["minima"]) == 2,
                             members=[dict(crest_index=item["crest_index"], relation=relation)
                                      for item, relation in group["members"]]))
    save_json(folder / "symmetry_groups.json", metadata)
    retained = [group["representative"] for group in groups]
    write_xyz(folder / "crest_pruned.xyz", retained)
    print(f"Symmetry pruning: {len(structures)} CREST structures -> {len(retained)} Gaussian inputs; "
          f"{sum(group['enantiomers'] for group in metadata)} mirror pairs.", flush=True)
    return retained


def populations(energies, multiplicities=None):
    delta = (np.asarray(energies) - min(energies)) * HARTREE_EV
    weights = np.exp(-delta / (lx.parser.BOLTZ_EV * TEMPERATURE))
    if multiplicities is not None:
        weights *= np.asarray(multiplicities)
    return delta, 100 * weights / weights.sum()


def conformer_sort_key(item):
    """One ordering for CSV groups and XYZ blocks, including energy ties."""
    return item["energy"], Path(item["source"]).name


def write_report(results, filename):
    results = sorted(results, key=conformer_sort_key)
    # Sum the weights of distinct handed minima, retaining separately calculated
    # Gibbs energies when both partners have Gaussian logs.
    def grouped_populations(key):
        reference = min(minimum[key] for item in results for minimum in item.get("minima", [item]))
        weights = [sum(np.exp(-(minimum[key] - reference) * HARTREE_EV /
                              (lx.parser.BOLTZ_EV * TEMPERATURE))
                       for minimum in item.get("minima", [item])) for item in results]
        delta = (np.array([item[key] for item in results]) - reference) * HARTREE_EV
        return delta, 100 * np.array(weights) / sum(weights)
    delta_e, pop_e = grouped_populations("energy")
    delta_g, pop_g = grouped_populations("gibbs")
    with open(filename, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Group", "E_Hartree", "DeltaE_eV", "PopE_300K_percent",
                         "G_Hartree", "DeltaG_eV", "PopG_300K_percent", "Gaussian_log",
                         "Multiplicity"])
        for index, item in enumerate(results):
            writer.writerow([index + 1, f"{item['energy']:.12f}", f"{delta_e[index]:.8f}",
                             f"{pop_e[index]:.5f}", f"{item['gibbs']:.12f}",
                             f"{delta_g[index]:.8f}", f"{pop_g[index]:.5f}",
                             Path(item["source"]).name, item.get("multiplicity", 1)])


def organize_outputs(folder):
    """Keep only the population report, manifest and unique ensemble at top level."""
    inputs = folder / "Inputs"
    cregen = folder / "CREGEN"
    inputs.mkdir(exist_ok=True)
    cregen.mkdir(exist_ok=True)
    for name in ("template.com", "search.json"):
        previous = folder / name
        if previous.is_file():
            previous.replace(inputs / name)
    for name in ("crest_reference.xyz", "crest_reoptimized.xyz",
                 "crest_reoptimized.xyz.sorted", "cregen.log",
                 "rejected_conformers.csv", "conformation.lx"):
        previous = folder / name
        if previous.is_file():
            previous.replace(cregen / name)
    return cregen


def classify_only(folder=".", crest="crest", rthr=0.125, ethr=0.05, charge=0, uhf=0, merge_mirrors=True):
    """Re-sort completed Gaussian opt/freq logs without repeating calculations."""
    folder = Path(folder).resolve()
    cregen_folder = organize_outputs(folder)
    files = sorted(set(folder.glob("Geometry-*.log")) | set((folder / "Geometries").glob("Geometry-*.log")))
    expected = list(folder.glob("Geometry-*.com")) + list((folder / "Geometries").glob("Geometry-*.com"))
    files = sorted(set(files) | {filename.with_suffix(".log") for filename in expected})
    results, rejected = [], []
    for filename in files:
        try:
            results.append(gaussian_result(filename))
        except (OSError, ValueError) as error:
            rejected.append((str(filename), str(error)))
    with open(cregen_folder / "rejected_conformers.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Gaussian_log", "reason"])
        writer.writerows(rejected)
    if not results:
        raise ValueError("No completed, verified minima. See CREGEN/rejected_conformers.csv.")
    if any(item["atoms"] != results[0]["atoms"] for item in results):
        raise ValueError("Gaussian outputs have inconsistent element sequences.")
    all_results = results
    groups = prune_symmetry(results, rthr, ethr, merge_mirrors)
    results = []
    membership = {}
    for group in groups:
        item = dict(group["representative"])
        minima = group["minima"][:]
        if len(minima) == 1 and merge_mirrors:
            # In an achiral equilibrium model every chiral minimum has an
            # equal-energy mirror minimum, even if sampling missed its partner.
            # Multiplicity must not depend on whether it happened to be sampled.
            reflected = dict(item, geometry=item["geometry"] * np.array([-1., 1., 1.]))
            if symmetry_rmsd(item, reflected) > rthr:
                minima.append(item)
        item.update(minima=minima, multiplicity=len(minima), enantiomers=len(minima) == 2,
                    partner_log=Path(minima[1]["source"]).name if len(minima) == 2 and
                    minima[1]["source"] != item["source"] else
                    ("inferred mirror partner" if len(minima) == 2 else ""))
        for member, relation in group["members"]:
            membership[member["source"]] = (item["source"], relation)
        results.append(item)
    energy_span = (max(item["energy"] for item in results) - min(item["energy"] for item in results)) * 627.509474
    write_xyz(cregen_folder / "crest_reference.xyz", [results[0]])
    write_xyz(cregen_folder / "crest_reoptimized.xyz", results)
    with TemporaryDirectory(prefix="cregen-", dir=cregen_folder) as temporary:
        sort_folder = Path(temporary)
        write_xyz(sort_folder / "reference.xyz", [results[0]])
        write_xyz(sort_folder / "reoptimized.xyz", results)
        command = [crest, "reference.xyz", "--cregen", "reoptimized.xyz", "--ewin", str(max(6.0, energy_span + 1)),
                   "--rthr", str(rthr), "--ethr", str(ethr), "--temp", "300",
                   "--chrg", str(charge), "--uhf", str(uhf)]
        try:
            run_command(command, sort_folder, "cregen.log")
        finally:
            logfile = sort_folder / "cregen.log"
            if logfile.exists():
                (cregen_folder / "cregen.log").write_text(logfile.read_text(encoding="utf-8"), encoding="utf-8")
        # Standalone CREGEN 3.0.2 writes <input>.sorted.
        output = sort_folder / "reoptimized.xyz.sorted"
        (cregen_folder / "crest_reoptimized.xyz.sorted").write_text(output.read_text(encoding="utf-8"), encoding="utf-8")
        unique = read_xyz(output)
        representatives = []
        for structure in unique:
            energy = float(structure["comment"].split()[0])
            candidates = [item for item in results if item["atoms"] == structure["atoms"]
                          and abs(item["energy"] - energy) <= 1e-6]
            matches = [(aligned_rmsd(item["geometry"], structure["geometry"]), item)
                       for item in candidates]
            if not matches:
                raise ValueError(
                    f"CREGEN output has no source with matching elements and energy "
                    f"({energy:.12f} Hartree). Compare crest_reoptimized.xyz "
                    f"and crest_reoptimized.xyz.sorted.")
            rmsd, match = min(matches, key=lambda pair: pair[0])
            if rmsd > 5e-4:
                raise ValueError(
                    f"CREGEN output does not match a source after alignment "
                    f"(best RMSD {rmsd:.6f} Angstrom). Compare atom ordering/units in "
                    f"crest_reoptimized.xyz and crest_reoptimized.xyz.sorted.")
            representatives.append(match)
        representatives.sort(key=conformer_sort_key)
        for group, item in enumerate(representatives, 1):
            item["comment"] = (f"{item['energy']:.12f} Group={group} "
                               f"Gaussian_log={Path(item['source']).name}")
        write_xyz(folder / "conformers_unique.xyz", representatives)
    write_report(representatives, folder / "conformation.csv")
    retained = {item["source"] for item in representatives}
    group_details = {item["source"]: (index, item)
                     for index, item in enumerate(representatives, 1)}
    with open(folder / "conformers_manifest.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Gaussian_log", "E_Hartree", "G_Hartree", "status", "Group",
                         "Multiplicity", "Enantiomers", "Partner_Gaussian_log"])
        for item in all_results:
            representative, relation = membership[item["source"]]
            group, details = group_details.get(representative, ("", {}))
            writer.writerow([Path(item["source"]).name, item["energy"], item["gibbs"],
                             relation + ": " + Path(representative).name if representative in retained else
                             "removed by CREGEN (duplicate/rotamer/topology)", group,
                             details.get("multiplicity", ""),
                             ("yes" if details.get("enantiomers") else "no") if details else "",
                             details.get("partner_log", "")])
    print(f"{len(all_results)} verified minima -> {len(representatives)} conformer groups; {len(rejected)} rejected logs.", flush=True)
    return representatives


def run_workflow(gaussian_input, crest_batch, gaussian_batch, gaussian="g16", max_jobs=10,
                 workdir="Conformational", xtb="xtb", crest="crest", solvent=None,
                 rthr=0.125, ethr=0.05, merge_mirrors=True):
    template = lx.parser.read_gaussian_input(gaussian_input)
    template["gaussian"] = gaussian
    if gaussian not in ("g16", "g09") or max_jobs < 1:
        raise ValueError("Use g16/g09 and a positive maximum number of simultaneous jobs.")
    gaussian_routes(template)
    batches = [Path(crest_batch).resolve(), Path(gaussian_batch).resolve()]
    for batch in batches:
        if not batch.is_file() or "#SBATCH" not in batch.read_text(encoding="utf-8"):
            raise ValueError(f"Expected a SLURM batch script that executes bash \"$1\": {batch}")
    folder = Path(workdir).resolve()
    folder.mkdir(exist_ok=True)
    organize_outputs(folder)
    xtb_folder, crest_folder, gaussian_folder = [folder / name for name in ("xTB", "CREST", "Geometries")]
    for stage in (xtb_folder, crest_folder, gaussian_folder):
        stage.mkdir(exist_ok=True)
    saved_input = folder / "Inputs" / "template.com"
    input_text = Path(gaussian_input).read_text(encoding="utf-8")
    if saved_input.exists() and saved_input.read_text(encoding="utf-8") != input_text:
        raise ValueError("Existing search uses a different Gaussian input; choose a new workdir.")
    settings_file = folder / "Inputs" / "search.json"
    if settings_file.exists() and json.loads(settings_file.read_text()).get("solvent") != solvent:
        raise ValueError("Existing search uses a different xTB/CREST solvent; choose a new workdir.")
    for obsolete in (folder / "Inputs" / "Interrupted", folder / "Geometries" / "Interrupted"):
        remove_output(obsolete)
    if not saved_input.exists():
        temporary = saved_input.with_name("template.com.tmp")
        temporary.write_text(input_text, encoding="utf-8")
        temporary.replace(saved_input)
    save_json(settings_file, dict(charge=template["charge"], uhf=template["multiplicity"]-1,
                                 nproc=template["nproc"], solvent=solvent))
    solvation = ["--alpb", solvent] if solvent else []
    if sampling_finished(xtb_folder, "xtbopt.xyz", "xtb.log",
                         "normal termination of xtb", template["atoms"], single=True):
        print("Reusing completed xTB optimization.", flush=True)
    else:
        restart_stage(xtb_folder)
        print("Optimizing the starting geometry with GFN2-xTB locally...", flush=True)
        write_xyz(xtb_folder / "start.xyz", [dict(atoms=template["atoms"], geometry=template["geometry"])])
        xtb_command = [xtb, "start.xyz", "--gfn", "2", "--opt", "tight", "--chrg", str(template["charge"]),
                       "--uhf", str(template["multiplicity"]-1), "--parallel", str(template["nproc"])]
        run_command(xtb_command + solvation, xtb_folder, "xtb.log")
        if not sampling_finished(xtb_folder, "xtbopt.xyz", "xtb.log",
                                 "normal termination of xtb", template["atoms"], single=True):
            raise ValueError("xTB optimization incomplete; inspect xTB/xtb.log.")
    optimized = read_xyz(xtb_folder / "xtbopt.xyz")
    if sampling_finished(crest_folder, "crest_clustered.xyz", "crest.out",
                         "CREST terminated normally", template["atoms"]):
        print("Reusing completed CREST search.", flush=True)
    else:
        restart_stage(crest_folder)
        # A new CREST ensemble may assign entirely different conformer numbers.
        restart_stage(gaussian_folder)
        write_xyz(crest_folder / "start.xyz", optimized)
        crest_command = [crest, "start.xyz", "--gfn2", "-T", str(template["nproc"]), "--cluster",
                         "--chrg", str(template["charge"]), "--uhf", str(template["multiplicity"]-1)] + solvation
        script = crest_folder / "run_crest.sh"
        script.write_text("#!/bin/bash\nset -e\n" + " ".join(shlex.quote(arg) for arg in crest_command)
                          + " > crest.out 2>&1\n", encoding="utf-8")
        print("Submitting CREST to SLURM...", flush=True)
        run_slurm(batches[0], script, crest_folder, template["nproc"], template["mem"])
        if not sampling_finished(crest_folder, "crest_clustered.xyz", "crest.out",
                                 "CREST terminated normally", template["atoms"]):
            raise ValueError("CREST search incomplete; inspect CREST/crest.out.")
    structures = read_xyz(crest_folder / "crest_clustered.xyz")
    structures = sampling_prune(structures, crest_folder, rthr, ethr, merge_mirrors)
    files = make_gaussian_inputs(template, structures, gaussian_folder)
    statefile = gaussian_folder / "frequency_retries.json"
    retries = json.loads(statefile.read_text()) if statefile.exists() else {}
    pending = [name for name in files if name not in retries
               and not gaussian_finished(gaussian_folder / Path(name).with_suffix(".log"))]
    if pending:
        print(f"Submitting {len(pending)} unfinished Gaussian opt/freq jobs; "
              f"maximum {max_jobs} simultaneous jobs...", flush=True)
        failed = run_gaussian_jobs(template, batches[1], pending, gaussian_folder, max_jobs)
        if failed:
            print("Failed Gaussian jobs (excluded from the final ensemble): " + ", ".join(failed), flush=True)
    else:
        print("Initial Gaussian opt/freq jobs already finished.", flush=True)
    retry_frequency_checks(template, batches[1], files, gaussian_folder, max_jobs)
    signature = classification_signature(folder, files, crest, rthr, ethr, merge_mirrors)
    result = cached_classification(folder, signature)
    if result is None:
        result = classify_only(folder, crest, rthr, ethr, template["charge"], template["multiplicity"]-1, merge_mirrors)
        outputs = {name: file_digest(folder / name) for name in
                   ("conformation.csv", "conformers_manifest.csv", "conformers_unique.xyz")}
        save_json(folder / "CREGEN" / "classification.json", dict(signature=signature, outputs=outputs))
    else:
        print("Reusing completed classification and population reports.", flush=True)
    print(f"Search complete. Report: {folder / 'conformation.csv'}", flush=True)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gaussian_input", nargs="?")
    parser.add_argument("--batch", help="Shared SLURM script for CREST and Gaussian; auto-detected if omitted.")
    parser.add_argument("--crest-batch")
    parser.add_argument("--gaussian-batch")
    parser.add_argument("--gaussian", choices=("g09", "g16"), default="g16")
    parser.add_argument("--max-jobs", type=int, default=10)
    parser.add_argument("--workdir", default="Conformational")
    parser.add_argument("--xtb", default="xtb")
    parser.add_argument("--crest", default="crest")
    parser.add_argument("--solvent", help="Optional xTB/CREST ALPB solvent; Gaussian SCRF is preserved separately.")
    parser.add_argument("--rthr", type=float, default=0.125)
    parser.add_argument("--ethr", type=float, default=0.05)
    parser.add_argument("--keep-enantiomers", action="store_true",
                        help="Keep mirror partners separate (e.g. for chiral environments).")
    parser.add_argument("--classify-only", metavar="FOLDER")
    args = parser.parse_args(argv)
    try:
        if args.rthr <= 0 or args.ethr < 0:
            raise ValueError("RMSD threshold must be positive; energy threshold must be nonnegative.")
        if args.classify_only:
            config = Path(args.classify_only) / "Inputs" / "search.json"
            if not config.exists():
                config = Path(args.classify_only) / "search.json"
            settings = json.loads(config.read_text()) if config.exists() else {}
            classify_only(args.classify_only, args.crest, args.rthr, args.ethr,
                          settings.get("charge", 0), settings.get("uhf", 0), not args.keep_enantiomers)
        else:
            if not args.gaussian_input:
                parser.error("Provide a Gaussian input file.")
            shared_batch = args.batch
            if not (shared_batch or args.crest_batch or args.gaussian_batch):
                shared_batch = lx.tools.find_batch_script()
            crest_batch = args.crest_batch or shared_batch or args.gaussian_batch
            gaussian_batch = args.gaussian_batch or shared_batch or args.crest_batch
            run_workflow(args.gaussian_input, crest_batch, gaussian_batch, args.gaussian,
                         args.max_jobs, args.workdir, args.xtb, args.crest, args.solvent, args.rthr, args.ethr,
                         not args.keep_enantiomers)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Conformational search failed: {error}", file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
