# Conformational search with CREST and Gaussian

The supported LeoX workflows are conformational search and omega tuning. The
legacy spectrum, exciton, and distortion features remain available but are
deprecated. The previous temperature-ramped Wigner conformational search is
replaced by this workflow:

1. Read one Gaussian `.com`/`.gjf` input containing explicit Cartesian coordinates
   in Angstroms. Element symbols or atomic numbers are accepted. Extract its
   full route, charge/multiplicity, processor count, memory, and trailing input
   sections (for example Gen/GenECP basis sets and SCRF=Read data).
2. Optimize the starting structure locally with `xtb --gfn 2 --opt tight`.
3. Submit `crest start.xyz --gfn2 -T N --cluster` through SLURM and wait for it.
4. Read `crest_clustered.xyz`; generate one Gaussian optimization followed by a
   frequency Link1 job per structure. Each structure is submitted as a separate
   SLURM job, and the existing Watcher monitors both normal terminations.
5. Use standalone CREGEN on the verified optimized minima, then write unique
   geometries and their electronic-energy and Gibbs-energy populations at 300 K.

The original high-level method, basis, SCF/integration settings, solvent
settings, and optimization options are retained. Frequency options from the
template are used in the second job (default `freq=noraman`). The thermochemistry
temperature is set to 300 K. Each job has its own checkpoint. Frequency steps
read the optimized geometry/wavefunction from that checkpoint.

Supply a **single-job** Gaussian template; LeoX generates the second Link1 job.
Z-matrices, frozen-atom coordinate columns, fragment charge specifications,
checkpoint-only inputs, transition-state searches, and read-in optimization
constraints are rejected. Coordinates are not read from the basis-set section.
For numerical frequency calculations, the final optimized geometry is extracted
from the optimization step, rather than from subsequent displaced geometries.

## Setup

Install/update LeoX, and load xTB and CREST in the shell that launches the driver:

```bash
python -m pip install .
# Use your cluster's actual module names:
# module load xTB CREST
```

Provide two SLURM batch scripts. Both must set your partition/time/module/scratch
configuration and execute the supplied command script with:

```bash
bash "$1"
```

Examples are provided as `batch_examples/crest_batch.sh` and
`batch_examples/gaussian_conf_batch.sh`. Replace their partition and module
placeholders. These are the actual SLURM scripts, not the old `slurm.sh` wrapper.
Do not put another `sbatch` invocation inside them. LeoX handles submission.
Relative auxiliary files referenced by the Gaussian template (for example
`@basis.gbs`) must be accessible from the `Geometries` working directory; inline
basis/solvent data are copied automatically.

LeoX requests one shared-memory task with N CPUs, where N comes from `%nproc` or
`%nprocshared` (default 1). It also requests the Gaussian `%mem` value (default
1GB), converting Gaussian byte/word memory syntax to SLURM units. The scripts'
partition, time limit, and module configuration remain in effect. Choose these
resources to suit both programs. Gaussian `%mem` does not include all process
overhead, so adjust the template memory or site resource setup if needed.

## Run

Select option 6 from `lx`, or use the command line:

```bash
lx_conf_search molecule.com \
  --crest-batch crest_batch.sh \
  --gaussian-batch gaussian_conf_batch.sh \
  --gaussian g16 --max-jobs 10 --workdir Conformational
```

The command-line driver stays running while its SLURM jobs execute. Launch it
inside tmux/screen, or detach it explicitly:

```bash
nohup lx_conf_search molecule.com \
  --crest-batch crest_batch.sh \
  --gaussian-batch gaussian_conf_batch.sh \
  --workdir Conformational > conformational.log 2>&1 &
```

The interactive `tools.conformational()` path starts the driver in the background
and writes `<workdir>.log`. New searches require a new output directory to avoid
overwriting previous calculations. Each Gaussian job contains one opt/freq pair;
`--max-jobs` controls simultaneous jobs, not the number of conformers per job.

The Gaussian SCRF setting is preserved. xTB/CREST use gas phase by default; add
`--solvent toluene` (or another supported ALPB solvent name) to use ALPB in both
sampling stages. Gaussian and ALPB solvent names/models are not assumed to be
interchangeable. Charge and unpaired-electron count (multiplicity minus one) are
passed to xTB and CREST.

## Outputs

All files are placed below the chosen output directory:

| Path | Contents |
| --- | --- |
| `template.com`, `search.json` | Original template and sampling metadata |
| `xTB/` | Starting XYZ, xTB log, and optimized structure |
| `CREST/` | SLURM submission log, CREST output and ensembles |
| `Geometries/` | Gaussian inputs, logs, command scripts and SLURM completion markers |
| `crest_reoptimized.xyz` | All verified high-level minima, with electronic energies in Hartree |
| `conformers_unique.xyz` | Representatives retained after CREGEN deduplication |
| `cregen.log` | Final sorting output |
| `conformation.lx` | Electronic and Gibbs energies, relative energies, populations and source logs |
| `rejected_conformers.csv` | Failed/incomplete jobs, imaginary frequencies, or missing thermochemistry |
| `conformers_manifest.csv` | Verified minima and whether CREGEN retained each representative |

Final sorting uses `--cregen`, compatible with CREST 3.0.2. The energy window is
set wide enough to retain all verified minima before duplicate/topology filtering.
`--rthr` (default 0.125 Angstrom) and `--ethr` (default 0.05 kcal/mol) expose CREGEN
thresholds. CREGEN also uses rotational constants and topology; it is not a
pure symmetry-aware RMSD comparison. Structures related by equivalent-atom
permutations can survive its classical sorting. Inspect borderline structures
before changing thresholds. Final PCA/k-means representative sampling is not
applied, because that would discard distinct minima from the population report.

Both population columns use

\[
p_i = \frac{\exp[-(X_i-X_{\min})/(k_B\,300\,\mathrm{K})]}
{\sum_j \exp[-(X_j-X_{\min})/(k_B\,300\,\mathrm{K})]}.
\]

`PopE` uses electronic energies; `PopG` uses Gaussian's harmonic rigid-rotor
Gibbs energies computed at 300 K. The latter includes thermochemical corrections
but can be sensitive to low-frequency modes; no quasi-harmonic correction is
introduced. Collapsed starting structures are not counted as statistical
degeneracy. Populations are normalized over retained verified minima; excluded
jobs and the initial selection of CREST cluster representatives limit completeness.

## Re-sort without repeating calculations

After correcting failed jobs manually, or to try another sorting threshold:

```bash
lx_conf_search --classify-only Conformational --rthr 0.125 --ethr 0.05
```

This reads existing Gaussian opt/freq logs, retrieves charge/spin metadata from
`search.json`, and regenerates the report. It does not repeat xTB, CREST sampling,
or Gaussian calculations. SLURM errors, cancelled jobs and incomplete Gaussian
logs are detected by the watcher and excluded rather than left waiting forever.
Failed Gaussian jobs are not retried automatically by this workflow.

## Validation

```bash
python -m unittest discover -s tests -v
```

Tests cover input/settings preservation, two-step Gaussian jobs, simulated
xTB/CREST/SLURM execution, duplicate collapse, failed jobs, numerical-frequency
geometry extraction, thermochemical energies, memory conversion and populations.
They do not replace a real run with the site-installed chemistry programs.
