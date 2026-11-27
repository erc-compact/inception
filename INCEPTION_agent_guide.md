# INCEPTION guide for AI agents

INCEPTION version 1.0.15. Source: <https://github.com/erc-compact/inception>.

This file is for an AI agent that helps a user run INCEPTION. Read all of it before acting. It covers what
INCEPTION does, how to run each part, what to check, and what to ask the user. Commands assume the
repository root as the working directory. The guide for humans, with figures and background, is
`INCEPTION_human_guide.pdf` next to this file; point the user there for explanations.

## How to work with the user

Follow this order:

1. **Find the goal.** Ask which of the tasks below (*Tasks*) the user wants. Offer the list as options.
2. **Discover, do not ask.** Check the software and read the filterbank headers yourself (see *Discover
   before asking*). Tell the user what you found. Ask only for what you cannot find.
3. **Interview.** Ask the open questions for that task (see *Questions to ask*) in groups of two to four. Give a
   recommended default with every question and accept "default" as an answer. Do not choose scientific
   values (S/N, DM, period or orbit ranges, search parameters) on your own; propose them and let the user
   decide.
4. **Propose.** Write the plan file and the exact commands. Show them with the expected outputs, the disk
   space needed and a rough run time for every step, not only the injection (see *Cost of a run*). Wait
   for approval.
5. **Test small.** Run one filterbank or one injection first and verify it (see *Verifying an
   injection*). Show the result. Scale up only after the user agrees.
6. **Report.** Give the paths of every output. Keep the plan, the injection report and the commands
   together so the run can be reproduced.

Rules:

- Never change, move or delete the user's input data. INCEPTION only reads its input filterbank.
- Ask before you submit cluster jobs, start anything that runs longer than a few minutes, download files
  (ephemerides, container images), delete anything (including Nextflow's `work/`), or write outside the
  agreed output directory.
- Check free disk space before injecting. Every injection writes a full-size copy of the input filterbank.
- Whenever you wait for something that runs longer than a few minutes (an injection, a cluster job, a
  Nextflow run), set up a **heartbeat of at most one hour**: a scheduled wake-up, a timer or a background
  check with a time limit, whatever your environment offers. Never wait on a job without a time limit. At
  each heartbeat check that the run is progressing (*Pipeline set-up*, step 8). If it has not progressed
  since the last heartbeat, treat it as hung: find the cause and tell the user. Cancel a hung job only with
  the user's approval; for unattended runs, ask for that approval before you start.
- Do not edit INCEPTION's code to make a run work. Tell the user what failed and why.
- Check every plan against *Pitfalls* and *Known issues* before you run it.
- When you start from an existing plan, show the user its `DM_smear`, `DM_const` and `SNR_def` (or say
  they are unset, so the code defaults `off`, `exact` and `total` apply) and confirm them.
- If the user explicitly delegates the choices ("pulsars of your choice, start processing"), choose
  sensible values yourself, state them, and proceed. Still ask before deleting anything.
- Warn the user when a plan puts more than about 20 pulsars into one filterbank, and recommend fewer
  (about 10 or less) when the pulsars have similar DMs (a globular cluster), very high S/N or long periods
  (see *Pulsars per filterbank*).
- If the user runs into a problem you cannot solve (an error that looks like a bug, behaviour that
  contradicts this guide, a missing feature), suggest they contact the INCEPTION developer, Robert Senzel
  (rsenzel@mpifr-bonn.mpg.de), and help them collect the plan, the command, the log and the injection
  report to send.

## Tasks

| Task | What the user wants | Tool | Read |
|---|---|---|---|
| A | Add simulated pulsars to one or a few filterbanks | standalone injector | Standalone injector, Injection plans, Parameter reference |
| B | Clean test data: pulsars in Gaussian noise with a real header | standalone injector, `--generate_fb` | Standalone injector |
| C | Inject a pulsar from a TEMPO/TEMPO2 par file (known, new or made up; any PINT timing or binary model) | standalone injector, `pint_polycos` | PINT injections |
| D | Measure how well a search or folding setup recovers pulsars, with many injections | Nextflow: `INJECTOR`, `PEASOUP` or `PRESTO` | Nextflow pipelines |
| E | Fold known pulsars or candidates through a pipeline without injecting anything | plan with `SNR: 0` (follow-up mode) | Follow-up mode, Nextflow pipelines |
| F | Remove bright known pulsars from data before a search | Nextflow: Nullsar | Nullsar |

**Recommended set-up for injections** (tasks A, B, D): for anything beyond a quick single-file test,
recommend the `INJECTOR` pipeline with only one folding package (see *Turning steps on and off*):
PulsarX with `SNR_def: "total"`, or PRESTO with `SNR_def: "pulsed"` if PulsarX is not available. Advise
against DSPSR for S/N measurements; its S/N follows the injected S/N less closely.

Suggested opening question: "What would you like to do? (A) add simulated pulsars to a filterbank,
(B) make simulated test data in pure noise, (C) inject a pulsar from a par file, (D) test how well a search
pipeline finds pulsars, (E) fold known pulsars or candidates, (F) remove bright known pulsars from your
data."

## Discover before asking

**Repository.** `inception/` contains `injector/` (library, `SCRIPT_inject_pulsars.py`,
`example_plans/`), `search_pipeline_validator/` (Nextflow pipelines and their Python scripts) and
`nullsar/` (Nextflow pipeline and code). Run the standalone injector from the repository root.

**Python environment** (all tasks):

```bash
python3 -c "import numpy, scipy, astropy, sympy, pandas, matplotlib"
python3 -c "import pint; print(pint.__version__)"
python3 -c "import jplephem"
```

The first line is required. `pint-pulsar` is needed for tasks C and F and for `polycos`/`t2pred`
predictors. `jplephem` is needed to pass a JPL `.bsp` file to `--ephem`. Ask where a JPL ephemeris such as
`de440.bsp` is; if there is none, offer to download it from <https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/>
(ask first) or use `--ephem builtin`.

**Pipelines** (tasks D, E, F):

```bash
java -version
nextflow -version
singularity --version || apptainer --version
sinfo -s
```

Nextflow needs Java 17 or later; both can be installed in the user's home directory without root rights
(`curl -s https://get.nextflow.io | bash`, then move `nextflow` into `~/bin`). `sinfo` lists SLURM
partitions; for another batch system, `cluster.config` must be adapted. Ask for the paths of the container
images (`.sif`): Python with the INCEPTION dependencies, PulsarX, PRESTO, DSPSR/PSRCHIVE, PEASOUP and,
optionally, the PICS classifier. Only the images of the steps that run are needed.

**Internet access.** Ask whether the machine that runs the injector, and the compute nodes, can reach the
internet. This decides `--cache` and `--offline` (see *Caches, offline nodes, SLURM*).

**Filterbank header.** Run from the repository root:

```bash
python3 -c "
import sys
from injector.io_tools import FilterbankReader
fb = FilterbankReader(sys.argv[1], stats_samples=0)
print(fb.header)
print('n_samples', fb.n_samples, 'obs_len_s', fb.obs_len, 'ftop', fb.ftop, 'fbottom', fb.fbottom)
" /path/to/obs.fil
```

Report `telescope_id`, `nchans`, `nbits`, `tsamp`, the length, the band, `tstart` (MJD) and the pointing
(`src_raj`, `src_dej` are stored as numbers, for example `51406.7` means 05:14:06.7). Check:

- `nifs` is 1 (summed polarisations).
- `nbits` is 1, 2, 4, 8 or 16 (unsigned), or 8 or 16 with the header keyword `signed`.
- `nchans × nbits` is a multiple of 8 (each time sample fills whole bytes).
- `telescope_id` is supported: MeerKAT 64, Arecibo 1, Nançay 3, Parkes 4, Jodrell Bank 5, GBT 6, GMRT 7,
  Effelsberg 8, SRT 10, LOFAR 11, CHIME 20, FAST 21, MWA 30.
- For 2-bit data, read *2-bit data* and tell the user its caveats; 1-bit data are untested. Set
  `--gulp_size_GB` (or `gulp_size_GB` in the pipeline configs) to nbits/8 of its value: 0.0025 at 2 bits.
- Free disk space in the output directory is at least the file size (times the number of filterbanks
  written at once).

## Questions to ask

Ask only what is still open after *Discover before asking*. Recommended defaults are in brackets.

### All injection tasks (A to E)

1. **Data.** Which filterbank(s)? Are they raw or already RFI-cleaned? [must be raw, see *Pitfalls*]
   Coherently dedispersed, and at which DM? [`cDM` 0] Several files that form one observation and must be
   merged? (Merging is done by the pipelines; for the standalone injector, merge first.)
2. **Pulsars.** A fixed list, or random draws from distributions? How many per filterbank? [about 10, at
   most about 20; fewer for similar DMs, very high S/N or long periods] For each parameter, a value or a
   distribution: spin period or frequency, DM, S/N, duty cycle; binary orbit (orbital period in hours,
   companion mass or projected semi-major axis, eccentricity); line-of-sight acceleration or jerk; profile
   shape; spectral index [0]; scattering.
3. **Comparison software.** Which software will fold or search the data? This sets
   - `SNR_def`: [`total`] for PulsarX (S/N against an off-pulse baseline); `pulsed` for PRESTO (χ²) and
     Fourier searches. Advise against DSPSR for S/N comparisons: its S/N calibration is worse;
   - `DM_const`: [`historical`] matches TEMPO, TEMPO2, PINT, PRESTO and PulsarX. The code default is
     `exact`, which gives DMs about 0.01% off in those packages.
4. **Realism.** Intra-channel smearing [`DM_smear: "exact"`; the code default is `off`]. Position [beam
   centre from the header] or offsets from it with a beam FWHM.
5. **Reproducibility.** Fixed integer seed [yes] or a new draw every run (`"random"`).
6. **Outputs.** Output directory. Ephemeris type [`par`] or PulsarX candidate file. Phase predictor
   [`t2pred` for compact or strongly accelerated binaries, otherwise none].
7. **Solar System ephemeris.** [a JPL `de440.bsp` if available, otherwise `builtin`]

### Task B

- Noise mean and standard deviation: measure them from real data of the same setup (`fb_mean`, `fb_std`
  of `FilterbankReader(path, stats_samples=1e6)`, see *Discover before asking*) [`128,6` for 8-bit data only
  if nothing is known]. Tell the user which values you used; the report records them as `fb_mean`,
  `fb_sigma` and `generate_fb`. A noise level that differs from the real beams is easy to miss later.
- The input file only provides the header and length.
- Independent noise needs a different seed for every file: the same seed and number of processes give the
  same noise.

### Task C

The par file path. Check that it is valid at the epoch of the data. DM offset or extra spin or acceleration terms on top of the timing model? S/N,
profile and spectrum as for any pulsar.

### Task D

- Pipeline: `INJECTOR` (inject and fold only; recommended with only the PulsarX parfold), `PEASOUP` (GPU
  acceleration search) or `PRESTO` (`accelsearch`).
- Steps to run: filtool, parfold, candidate folding, classifier (see *Turning steps on and off*).
- Several epochs of the same pulsars? (see *Multi-epoch injections*)
- Target search S/N? It differs from the injected S/N (see *Pitfalls*); propose a pilot run.
- Number of injections (`start`..`end`) and how many at once [`batch_size` 25].
- The data list: filterbanks per entry, coherent DM per entry, random or cycling choice per injection.
- Search settings: DM plan per time-downsampling factor, segment plan (full observation, halves, ...),
  acceleration range (PEASOUP) or `zmax`/`wmax` (PRESTO), detection threshold, harmonics, FFT sizes. Full DM
  plan, or only the injected DMs (`inj_DM`, much faster)?
- Telescope settings: channel mask (MHz), birdies (Hz), PulsarX template.
- Cluster: partitions, GPU partition for PEASOUP, node-local scratch, directories to bind into the
  containers, image paths, account or QOS options.
- Which products to keep or delete; disk quota.

### Task F

- Filterbanks and their tags (files with the same tag form one observation).
- Par files of the pulsars to remove (accurate at the epoch of the data).
- S/N limit below which a pulsar is left alone [15].
- filtool options [`--td 4 -v`; downsampling is optional].
- Gaussian components per profile [2, or `"auto"` for complex profiles].
- PulsarX template, number of sub-integrations [64], cluster config.

## Standalone injector

### Injector command

```bash
mkdir -p out
python3 injector/SCRIPT_inject_pulsars.py --signal=plan.json --fb=obs.fil --output=out/ --ephem=de440.bsp --ncpu=8
```

| Option | Default | Meaning |
|---|---|---|
| `--signal` | required | Injection plan (JSON), see *Injection plans*. |
| `--fb` | required | SIGPROC filterbank. Only read. |
| `--output` | required | Output directory. Must exist; the injector does not create it. |
| `--ephem` | `builtin` | JPL `.bsp` file (needs `jplephem`), an Astropy ephemeris name (`de440`) or `builtin`. |
| `--ncpu` | 1 | Processes; the data are split in time. Use the number of cores. |
| `--gulp_size_GB` | 0.01 | Block size per read and write, in file bytes. Below 8 bits a block holds 8/nbits times more samples and memory: scale by nbits/8 (0.0025 at 2 bits). |
| `--stats_samples` | 1e6 | Samples used to measure channel means and standard deviations. |
| `--cache` | none | Directory for downloads (Astropy IERS tables, PINT clock files). |
| `--offline` | off | Never download. Uses cached or bundled tables and built-in telescope positions. |
| `--generate_fb` | none | `mean,std`: write Gaussian noise instead of copying the data, then inject (task B). |

Facts:

- All pulsars of a plan are summed and injected in one pass. Put pulsars meant for the same data in one
  plan; do not inject one after another.
- The output filterbank is a full copy of the input (same size).
- If every pulsar has `SNR: 0` and `--generate_fb` is not set, nothing is injected and the input is
  hard-linked (or copied, across file systems) to the output name. That output is the input file: never
  modify it in place.
- With `--generate_fb`, fresh noise is written even if every pulsar has `SNR: 0`. One pulsar with
  `SNR: 0` plus `--generate_fb` makes a pure-noise file.
- The script prints all pulsar parameters and their defaults when it starts.
- Relative paths inside a plan (`profile`, `gain_map`, `pint_polycos`) are resolved from the current
  directory. Use absolute paths when unsure.

### Injector outputs

| File | Content |
|---|---|
| `<fb stem>_<injection_id>.fil` | Injected filterbank, input header. |
| `report_<injection_id>_<seed>.json` | Observation metadata and every resolved pulsar parameter, including random draws. Can be used as a plan to repeat the injection (in-built timing model only). |
| `<ID>.par` | TEMPO par file of the injected pulsar (`create_parfile: "par"`). |
| `<ID>.candfile` | PulsarX candidate file (`create_parfile: "pulsarx"`). |
| `<ID>_predictor.t2pred`, `<ID>_predictor.polycos` | Phase predictor from the injected model (`create_predictor`). |
| `<ID>.polycos` | PINT polycos (task C). |

### Ephemeris files

How plan parameters map to the written ephemerides:

| Plan | Par file (`par`) | PulsarX candfile (`pulsarx`) |
|---|---|---|
| `RAJ`, `DECJ`, `PMRA`, `PMDEC`, `POSEPOCH` | same | not used (topocentric fold) |
| `DM` | `DM` | `DM` |
| `P0`/`F0`, `P1`/`F1`, ... | `F0`, `F1`, ... | topocentric `F0` at mid-observation; `F1` = 0 |
| `PEPOCH` | `PEPOCH` (moved to the acceleration epoch if `A0`... set) | mid-observation (PulsarX default) |
| `A0`, `A1`, ... | converted into `F1`, `F2`, ... | `accel` = `A0`, `F2` from `A1` |
| `binary_period` (h) | `PB` (d), `BINARY BT` | `Pb` (d) |
| `x` or masses + `inc` | `A1` | `A1` |
| `ecc` | `ECC` | `ECC` |
| `AoP` (deg) | `OM` (deg) | `OM` (radians, as PulsarX's candfile reader expects; `--OM` on its command line is degrees) |
| `T0` | `T0` | `T0` in the telescope frame |
| `fold_harmonic` | frequencies / harmonic | `F0` / harmonic |

- Par file also has `UNITS TDB`, `CLK TT(BIPM)`, `EPHEM DE421`, `TZRMJD` = `PEPOCH` at the barycentre and
  infinite frequency (top of band with `DM_ref: "top"`), a placeholder `PSR` name, no post-Keplerian terms.
- The par DM is the injected DM; use `DM_const: "historical"` so folding software reproduces the delays.
- The par file is barycentric: for `frame: "topo"` fold with the candfile or a predictor.
- The candfile has no spin-down and refers to mid-observation: keep `PEPOCH`, `ACCEPOCH` at defaults.
- PINT pulsars: the input par file is copied with the DM offset and extra spin terms added.
- Predictors (`create_predictor`) are fitted to the injected phase over the observation ± 60 s (target
  better than 1e-6 turns at band centre; the injector prints the reached accuracy). `t2pred`: PulsarX
  `--t2pred <file>`, DSPSR `-P <file>`. `polycos`: PRESTO `prepfold -polycos <file> -psr <name>`, DSPSR
  `-P <file>`. Use them for compact orbits, orbits plus accelerations, or PINT models with extra terms.

### 2-bit data

- Injecting into existing 2-bit data uses two approximations: the outer levels are treated like all others,
  although a digitiser puts the noise tails into them (injected pulsar slightly too bright), and the noise
  is measured from the stored values, which look narrower than the real noise (slightly too faint).
- They cancel for a well-set digitiser: three thresholds at the noise mean and about one standard deviation
  above and below it. Tested on 2-bit Parkes data: injected pulsars fold at the same S/N, to within one per
  cent, as pulsars digitised together with the noise.
- If the thresholds lie much closer to or further from the mean, the injected S/N can be off by a few per
  cent. The share of samples at each level, per channel, shows where the thresholds lie.
- `--generate_fb` is not affected: the pulsar is added before rounding.
- 2-bit data lose S/N to the rounding compared with 8-bit, for injected and real pulsars alike; the injected
  S/N refers to the signal before digitisation, so injections measure this loss directly.

### Caches, offline nodes, SLURM

- The injector downloads Astropy IERS tables (when older than 30 days) and, for PINT, observatory clock
  files.
- Nodes without internet: run once with `--cache /shared/dir` where there is internet, then use
  `--cache /shared/dir --offline` on the nodes.
- `--offline` without a cache uses Astropy's bundled IERS tables (fine for searches). PINT injections then
  fail at once because the clock files are missing.
- `--ephem builtin --offline` runs reliably on compute nodes without internet for pulsars with the in-built
  timing model (no PINT).
- `injector/SLURM_inject_pulsars.sh` is an example batch script (copies the data to node-local scratch,
  runs the injector in a Singularity image, copies the results back).
- Measured cost on an HPC system (7.6 GB 8-bit beam, `--generate_fb`, 24 processes): pure noise about
  5 min, 10 pulsars 16–19 min (5–6 CPU-h), peak memory 4.5–9 GB. 24 cores and about 16 GB were enough;
  whole nodes are rarely needed. At 2 bits the default gulp made 24 processes exceed 64 GB and swap; with
  `--gulp_size_GB=0.0025` a 2.5 GB 2-bit file took 6–7 min and 27–35 GB.

## Injection plans

### Plan structure

```json
{
    "psr_global": {"injection_id": "my_run", "global_seed": 1234, "DM_const": "historical", "DM_smear": "exact"},
    "pulsars": [
        {"ID": "psr_A", "P0": 0.004, "DM": 80, "SNR": 20},
        {"ID": "psr_B", "F0": 333.3, "F1": -1e-15, "DM": 81.5, "SNR": 12, "binary_period": 5, "M2": 0.2}
    ]
}
```

- `injection_id` names the output filterbank and report [`inj_<seed>`].
- `global_seed`: an integer (reproducible) or `"random"`.
- Any other key in `psr_global` is a default for every pulsar; a pulsar's own key overrides it.
- Each pulsar needs an `ID`, `P0` or `F0` (unless `pint_polycos` is given) and `SNR`.

### Random draws

| Draw | Syntax |
|---|---|
| uniform | `{"rng": "uniform", "low": 1, "high": 30}` |
| log-uniform | `{"rng": "loguniform", "low": 0.001, "high": 0.03}` |
| normal | `{"rng": "normal", "mean": 0.3, "sigma": 0.05}` |
| choice | `{"rng": "choice", "samples": [0, 0.1, 0.5], "weights": [7, 2, 1]}` (weights optional) |

- Optional `"units"`: a number multiplies `low`, `high`, `mean`, `sigma`; `"T_obs"` uses the observation
  length in hours (for `binary_period`); `"dt"` the sampling time.
- `"ID": {"replicate": 5, "prefix": "msp"}` makes pulsars `msp_R0` to `msp_R4`, each with its own seed.
- **Every pulsar and every replicate gets its own seed**, drawn from the global seed in plan order, so
  their draws differ. A `seed` given in a pulsar entry overrides the drawn one.
- **A draw inside `psr_global` is made once per injection and shared by all its pulsars.** For example, a
  DM draw in `psr_global` gives every pulsar the same DM.
- Keep the number of pulsars per filterbank small (see *Pulsars per filterbank*). Build larger
  populations with many injections in the pipelines (see *Nextflow pipelines*).

Example population, ten pulsars (`injector/example_plans/07_population.json`). Use it only as a starting
point to propose to the user:

```json
{"ID": {"replicate": 10, "prefix": "msp"},
 "P0": {"rng": "loguniform", "low": 0.0015, "high": 0.03},
 "DM": {"rng": "uniform", "low": 0, "high": 300},
 "SNR": {"rng": "uniform", "low": 8, "high": 20},
 "duty_cycle": {"rng": "loguniform", "low": 0.02, "high": 0.3},
 "binary_period": {"rng": "uniform", "low": 0.5, "high": 10, "units": "T_obs"},
 "M2": {"rng": "normal", "mean": 0.3, "sigma": 0.08},
 "ecc": {"rng": "choice", "samples": [0.0, 0.1, 0.5], "weights": [0.7, 0.2, 0.1]}}
```

### Pulsars per filterbank

Every injected pulsar adds signal to the same data, and they are all searched together. Too many in one
file change the data's statistics and make the pulsars compete with each other in the candidate lists, so
the results no longer describe a single pulsar in that data.

- Default: about 10 pulsars per filterbank. Never more than about 20; warn the user if a plan asks for
  more, including through `replicate`.
- Fewer still when the pulsars
  - have similar DMs, as in a globular cluster: their candidates fall in the same DM trials and their
    harmonics can be sifted together;
  - have very high S/N: each one then changes the data and fills the candidate lists more;
  - have long periods: each one occupies more of the low-frequency spectrum with its many harmonics and has
    long pulses that affect running-baseline and normalisation steps.
- For a large population, run the same plan on many filterbanks or many injection numbers with the
  pipelines; each injection draws a new set of pulsars.

### Example plans

`injector/example_plans/` (run from the repository root; each was tested):

| Plan | Shows |
|---|---|
| `01_minimal.json` | Smallest valid plan. |
| `02_binaries.json` | Orbits from masses and inclination, from companion mass, from `x`; eclipse light curve; `t2pred`. |
| `03_double_pulsar.json` | Pulsar B as the companion of A (`double_pulsar`). |
| `04_cluster_accelerations.json` | Cluster pulsars, per-pulsar DM drawn around the cluster DM, `A0`, `A1`, `presto_z`. |
| `05_custom_profiles.json` | Fixed, randomised and file-based profiles. |
| `06_propagation.json` | Scattering, DM–scattering relation, coherently dedispersed data (`cDM`). |
| `07_population.json` | Replicated random population. |
| `08_known_pulsar_pint.json` | PINT injection, with a DM offset and extra acceleration. Replace the par path. |
| `09_followup.json` | Follow-up mode, `SNR: 0`: a known pulsar (par file) and an accelerated candidate. Replace the par path. |

## Parameter reference

"Fraction" means an MJD or a number between 0 and 1, read as a fraction of the observation after its start.

| Parameter | Default | Meaning |
|---|---|---|
| **Identity and frame** | | |
| `ID` | required | Unique name, used in file names. May be a replicate block. |
| `seed` | derived | Seed for this pulsar's draws. |
| `frame` | `bary` | `bary`: real pulsar, barycentred. `topo`: periodic at the telescope. |
| **Position** | | |
| `RAJ`, `DECJ` | header | hh:mm:ss, dd:mm:ss. Default: `src_raj`, `src_dej` of the filterbank. |
| `separation`, `position_angle` | 0 | Offset from the beam centre (arcmin) and its angle (deg, east of north). |
| `beam_fwhm` | 0 | Beam FWHM (arcmin); if set, the S/N is scaled by a Gaussian beam at the offset. |
| `PMRA`, `PMDEC` | 0 | Proper motion (mas/yr, PMRA includes cos Dec). |
| `DIST` / `parallax` | 10 kpc / none | Distance (kpc) or parallax (mas). |
| `RV` | 0 | Radial velocity (km/s, positive receding). |
| `POSEPOCH` | obs. start | MJD or fraction. |
| **Spin** | | |
| `P0` or `F0` | required | Period (s) or frequency (Hz). |
| `P1`, `P2`, ... / `F1`, `F2`, ... | 0 | Derivatives of any order; the two kinds can be mixed. |
| `PEPOCH` | obs. start | MJD or fraction. |
| `phase_offset` | 0 | Spin phase at PEPOCH. |
| **Acceleration** | | |
| `A0`, `A1`, `A2`, ... | none | Line-of-sight acceleration (m/s²), jerk (m/s³), ...; positive away from the observer. |
| `ACCEPOCH` | 0.5 | Reference epoch of the accelerations, as a fraction of the observation length after PEPOCH. |
| `presto_z`, `presto_w` | none | Acceleration and jerk in PRESTO units, instead of `A0`, `A1`. `z = -F0 A0 T_obs² / c`. |
| **Signal and emission** | | |
| `SNR` | required | S/N at emission, before propagation. `0`: do not inject (follow-up). Negative: subtract. |
| `SNR_def` | `total` | `total` or `pulsed` (*Questions to ask*, comparison software). |
| `duty_cycle` | 0.1 | FWHM of the default Gaussian profile, in periods. |
| `profile` | `default` | `default`; a dict of Gaussian components `{"phase": [...], "duty_cycle": [...], "amp": [...]}` (or `"n_components": n` with draws per key); a `.npy` file (1D, or 2D `nchan × nbin`); an EPN database `.txt` profile (<https://psrweb.jb.man.ac.uk/epndb/>). |
| `spectrum` | 0 | Spectral index relative to the centre frequency, or a `.npy` 1D spectrum of any length, sampled evenly from the bottom to the top band edge (**lowest frequency first**, whatever the file's channel order) and interpolated to the channels; non-negative values; only the shape matters, the S/N is kept. Differs from `gain_map` with `gain_axis: "freq"` (`nchan` values in the file's channel order, not normalised, changes the S/N). |
| `gain_map` | none | `.npy` gain map: 1D in time or frequency, or 2D `nchan × n_subints` (channel order of the file). Not normalised by the injector. |
| `gain_axis` | `time` | Axis of a 1D map: `time` or `freq` (then `nchan` values). |
| `micro_structure` | 0 | Experimental micro-structure timescale (µs). |
| **Propagation** | | |
| `DM` | 0 | pc cm⁻³. With `pint_polycos` it is an offset added to the par-file DM. |
| `cDM` | 0 | Coherent dedispersion DM of the data; smearing uses abs(DM − cDM). |
| `DM_smear` | `off` | `exact` (recommended), `approx` (default profile only) or `off`. Use `approx`/`off` only to test their effect. |
| `sample_smear` | `off` | `on`: average the profile over one sampling time, as the backend does (pulses narrower than a sample lose S/N). Matters only for pulses not much wider than a sample (off: 0.1-ms pulses at 1-ms sampling fold about 2.7x too strong, 40-µs at 0.1 ms 1.4x; wide pulses < 0.5%); keep `off` for profiles that already contain it (Nullsar's fold-based models). |
| `DM_ref` | `inf` | Reference frequency of the delays: `inf` or `top`. |
| `DM_const` | `exact` | `exact` (4148.806), `historical` (4149.378, TEMPO/TEMPO2/PINT/PRESTO/PulsarX) or a number. |
| `scattering_time` | 0 | ms at the centre frequency. |
| `scattering_index` | −4 | Frequency scaling of the scattering time. |
| `scattering_model` | none | `TPA` or `CORDES`: scattering time from the DM. |
| `scattering_model_sigma` | 0 | Offset from the `CORDES` relation, in units of its scatter. |
| **Binary orbit** (BT, no post-Keplerian terms) | | |
| `binary_period` | none | Orbital period in **hours**. Setting it turns the orbit on. |
| `x` | derived | Projected semi-major axis (light-seconds). |
| `M1`, `M2`, `inc` | 1.4, derived, 90 | Masses (solar) and inclination (deg). Give three of `x`, `M1`, `M2`, `inc`, or only `x`, or only `M2`. |
| `M_func` | 0 | Mass function (solar), instead of `x`. |
| `ecc`, `AoP` | 0, 0 | Eccentricity, argument of periastron ω (deg). The same angle as `OM` (longitude of periastron) in par files, and written there as `OM`. |
| `T0` | obs. start | Periastron epoch, MJD only. |
| `double_pulsar` | none | ID of another pulsar in the plan; this one becomes its companion (same orbit, masses swapped, AoP + 180°). |
| **Output files** | | |
| `create_parfile` | `par` | `par` or `pulsarx`. |
| `create_predictor` | `none` | `t2pred` or `polycos`. |
| `fold_harmonic` | none | Write the ephemeris at this multiple of the period. |
| **PINT** | | |
| `pint_polycos` | none | `.par` file (converted with PINT) or an existing `.polycos` file. Replaces the in-built timing model. |
| `pint_N`, `pint_T` | 12, 5 | Coefficients per polyco segment and segment length (minutes). Shorter or more for compact binaries. |
| `pint_nfreq` | 8 | Frequencies across the band for polycos, besides the centre; doubled automatically (up to 64) while the injected phase deviates from PINT by more than 10 ns (1: centre only, never increased). |

Notes:

- The injector fixes the S/N, not the flux: a broader pulse at the same S/N has a lower peak. Smearing,
  scattering, masked channels and RFI lower the recovered S/N, as for a real pulsar.
- Accelerations are written to the par file as spin derivatives (`F1 = -F0 A0 / c`, `F2 = -F0 A1 / c`).
  For compact orbits use the binary model, not large accelerations and jerks. Accelerations can be
  combined with an orbit or with a PINT ephemeris.
- Post-Keplerian effects need a par file through PINT (see *PINT injections*).

## PINT injections

### PINT plan

```json
{
    "psr_global": {"injection_id": "known", "global_seed": 1},
    "pulsars": [
        {"ID": "J1748-2446A", "pint_polycos": "/abs/path/J1748-2446A.par", "SNR": 40, "duty_cycle": 0.08, "spectrum": -1.6}
    ]
}
```

The command is the same as for any plan. Add `--cache` (PINT downloads clock files). PINT uses the par
file's `EPHEM`; pass the same kernel with `--ephem` (INCEPTION copies it into the astropy cache, so it works
offline; otherwise PINT downloads it). If PINT
cannot load it, the injector falls back to the `--ephem` kernel (written to `out/<ID>.par`) or to the
built-in ephemeris (about 10 µs, not timing-grade), with a WARNING.

### Requirements and checks

- `pint-pulsar` installed; clock files reachable or cached.
- The par file is valid at the epoch of the data.
- Several epochs: the par file needs `TZRMJD` (see *Multi-epoch injections*).
- PINT injections are timing-grade from version 1.0.14 (validated end to end with DSPSR, PSRCHIVE, PINT
  and TEMPO2). If the user times the injected data, give them the settings in *Timing injected data*.

### What happens and what is written

1. PINT writes polycos for the telescope at the band centre and at `pint_nfreq` (8) further frequencies
   across the band, covering the observation ±10 minutes, segments of `pint_T` minutes with `pint_N`
   coefficients, with the par file's `EPHEM`, `CLK` and `PLANET_SHAPIRO`. Each set is checked against PINT;
   segments are halved until the error is below 1e-6 turns (also for compact and relativistic binaries).
   Between the frequencies the delay is
   interpolated (ν⁻² and ln ν, ln² ν terms fitted, spline in ln ν for the rest), so every channel follows
   the full model (Earth's Doppler on the DM delay, DMX, FD, solar wind) for any channel count and
   bandwidth. `pint_nfreq: 1`: centre only, other channels follow the par-file DM alone.
   Then the injected phase is compared with PINT in 64 channels spread over the band (all if fewer) × 9
   times; above 10 ns `pint_nfreq` is doubled (up to 64, only while that halves the deviation). Log line:
   `PINT polycos for <ID> (<EPHEM>, <CLK>): ... Injected phase against PINT over <n> of <N> channels:
   largest <x> ns, rms <y> ns`; WARNING if 10 ns is not reached. Tested at DM 200–300 from LOFAR HBA to
   the Parkes UWL, 64–4096 channels: at most 0.9 ns, 0.2 ns rms with the default.
2. The pulsar is injected with them; they are saved as `out/<ID>.polycos`.
3. `out/<ID>.par` is written: the input par with the plan's DM offset and extra terms. Fold with it.

- `DM` is an offset from the par-file DM [0].
- `F0`, `F1`, ... and `A0`, `A1`, ... are added to the timing model.
- `create_predictor: "t2pred"` writes a TEMPO2 predictor of the final model (PulsarX folds it directly).
- `pint_polycos` can point at an earlier run's `.polycos` file to skip PINT (no par file is written then).

## Nextflow pipelines (tasks D and E)

### The three pipelines

Folder: `search_pipeline_validator/NEXTFLOW_pipelines/pipeline_general_<NAME>/` with `main.nf`,
`processes.nf`, `nextflow.config`, `cluster.config` and `configs/pipeline_configs.json`.

| NAME | Steps in `main.nf` |
|---|---|
| `INJECTOR` | `injection` → `pulsarx_parfold`, `presto_parfold`, `dspsr_parfold` (in parallel) → `collector` |
| `PEASOUP` | `injection` → `filtool` → `peasoup_setup` → `peasoup_search` (one job per DDplan entry × segment) → `match_candidates` → `pulsarx_candfold` (one per matched segment) → `classifier` → `collector`; `pulsarx_parfold` as a side branch |
| `PRESTO` | `injection` → `rfifind` → `presto_ddplan_setup` → `presto_dedisperse` (per DDplan entry) → `presto_search_setup` → `presto_fft`, `presto_accelsearch` (per DDplan entry × segment) → `presto_cleanup` → `presto_sift` → `match_candidates` → `presto_candfold` → `classifier` → `collector`; `presto_parfold` as a side branch |

How they work:

- Each injection number (`start`..`end`) is one injection. Channels carry only the injection number (after
  a fan-out, the number and a tag); data are written to `output_dir`. Each process emits the value it
  received once its job is done.
- Fan-out: a setup process writes a plan file, one line per job (for example
  `inj_000001_DDPLAN_2_SEG_1_2`). `expand_plan` makes one job per line; `collapse_tag` waits for all of
  them and passes the injection number on.
- The pipeline replaces the plan's `global_seed` with a seed derived from the injection number, and sets
  `cDM` from the data entry. The same `start`..`end` with the same plan and data list reproduces the same
  injections.
- **Injection numbers must be unique across a whole campaign**, not only within one run. The number alone
  sets the seed, so reusing it in another run, output folder or plan repeats the same random draws and,
  with `generate_fb`, the same noise. Keep a table of the number ranges used (for example stage A
  101–200, stage B 201–300).
- A run takes one plan. For several plans, launch one run per plan with its own number range, and submit
  each head job to the batch system instead of starting many Nextflow heads on a login node.
- Process labels for `cluster.config` are in `processes.nf`. PRESTO's two setup processes share the label
  `presto_setup`.

### Pipeline set-up

1. Copy the pipeline folder to the run location (scratch or project space). Nextflow writes `work/` where
   it is launched.

   ```bash
   cp -r inception/search_pipeline_validator/NEXTFLOW_pipelines/pipeline_general_PEASOUP /scratch/$USER/run1
   ```

2. Edit `nextflow.config` (`params`):

   | Key | Meaning |
   |---|---|
   | `injection_dir` | Folder with the plan and the pipeline config. |
   | `output_dir` | Results, one `inj_XXXXXX/` per injection. |
   | `pipeline_code` | Absolute path to `inception/search_pipeline_validator`. |
   | `tmp_dir` | Node-local scratch for every job. |
   | `injection_plan`, `config_params` | `${params.injection_dir}/injection_plan.json`, `.../pipeline_configs.json`. |
   | `*_image` | Container images of the steps that run. |
   | `start`, `end` | Injection numbers. |
   | `batch_size` | Maximum injections running at once (limits injected filterbanks on disk). |

3. Edit `cluster.config`: a `singularity` (or `apptainer`) block with `enabled = true`, `autoMounts = true`
   and `runOptions = '--nv -B /scratch -B /data ...'` (every path used must be bound; `--nv` for the PEASOUP
   GPU), and a `process` block with `executor = 'slurm'`, `queue`, `errorStrategy = 'retry'`, `maxRetries`,
   and one `withLabel: '<label>' { cpus, memory, time, queue, clusterOptions }` block per label. Useful
   patterns: `memory = { 36.GB * task.attempt }` to grow on retry; `errorStrategy = 'ignore'` for optional
   steps such as the classifier; `executor { queueSize = 200 }`. Start from the file in the pipeline folder
   and ask the user for the partitions and limits. Rough starting resources: injection 24 CPUs and 16 GB
   (measured for 7.6 GB beams, see *Caches, offline nodes, SLURM*), filtool 2 CPUs, PulsarX folds 2–10
   CPUs, PEASOUP search 1 CPU + 1 GPU, setup, matcher and collector 1 CPU.

4. Put the injection plan at `configs/injection_plan.json`. Use absolute paths inside it.

5. Edit `configs/pipeline_configs.json` (next section), in this order: `data`, `injection_args`, one section
   per step in `main.nf`, `collector_args`.

6. Test with `start = 1`, `end = 1`. Check the output folder and the folds (see *Verifying an injection*).

7. Launch from the pipeline folder. The head job must outlive your session; submit it as a job:

   ```bash
   sbatch --time=7-00:00:00 --cpus-per-task=1 --mem=8G --wrap="nextflow run main.nf -c cluster.config"
   ```

   (or run `nextflow run main.nf -c cluster.config` inside `tmux` or `screen`). Add `-resume` to continue
   an interrupted run without repeating finished jobs.

8. Monitor with a heartbeat of at most one hour (see *How to work with the user*). At each heartbeat check:
   - `squeue -u $USER`: the head job is running; states and pending reasons of the other jobs (`Priority`,
     `Resources` or a per-user job limit are normal queuing, not a hang);
   - the head job's output and when `.nextflow.log` last changed (Nextflow writes to it whenever a job is
     submitted or finishes);
   - the number of finished injections, for example `ls <output_dir>/inj_*/report_*.json | wc -l`, compared
     with the previous heartbeat.

   A job running more than about three times as long as the same step in the test, or no change in
   `.nextflow.log` for an hour while jobs are running, means a hang. Read that job's work directory, tell
   the user, and with their approval cancel it (`scancel`); the templates' `errorStrategy = 'retry'` then
   reruns it, otherwise continue with `-resume`. Set each process's `time` in `cluster.config` to about
   three times its test run time (not less than
   30 minutes), so that SLURM also ends hung jobs by itself.

   A failed job's error names its work directory; read `.command.sh` (exact command), `.command.err` and
   `.command.log` there.

9. After a successful run, with the user's approval, delete `work/` and `.nextflow.log*`. Keep them while
   `-resume` may still be needed.

### `pipeline_configs.json`

The `configs/pipeline_configs.json` in each pipeline folder is a complete, working example for MeerKAT
L-band data (paths written as `...`). Start from it, replace the paths and ask the user for the
telescope-specific values: channel masks (MHz), PulsarX template, birdies, `ddplan`, `segment_plan`,
`filtool` options, `dspsr` telescope (`-k`). Notes on the shipped templates:

- INJECTOR folds with all three packages; for the recommended PulsarX-only set-up, comment out the other
  folds and set `presto_parfold`, `dspsr_parfold` to `false` in `collector_args`.
- PEASOUP: 3 downsampling factors × (1 + 2) segments = 9 search jobs per injection; a segment may set
  `ram_limit_gb` (default 5 × file size); `collector_args.classifier` is `false` although `main.nf` runs the
  classifier.
- PRESTO: `inj_DM` is `true` (only injected DMs searched); `filtool_args` only matters with
  `rfi_cleaner: "filtool"` or `fold_mode: "filtool"`, which also need the `filtool` step called in
  `main.nf`.

Sections present per pipeline (read the template in the folder for every key and its current value):

| Section | Step | Pipelines |
|---|---|---|
| `data` | injection | all |
| `injection_args` | injection | all |
| `pulsarx_parfold_args` | PulsarX parfold | INJECTOR, PEASOUP |
| `presto_parfold_args` | prepfold parfold | INJECTOR, PRESTO |
| `dspsr_parfold_args` | DSPSR parfold | INJECTOR |
| `filtool_args` | filtool | PEASOUP, PRESTO |
| `peasoup_args` | PEASOUP setup and search | PEASOUP |
| `presto_rfi_args` | rfifind | PRESTO |
| `presto_search_args` | DDplan, dedispersion, FFT, accelsearch | PRESTO |
| `presto_sift_args` | sifting | PRESTO |
| `candidate_matcher_args` | matching | PEASOUP, PRESTO |
| `pulsarx_candfold_args` | candidate folding | PEASOUP |
| `presto_candfold_args` | candidate folding | PRESTO |
| `classifier_args` | PICS | PEASOUP, PRESTO |
| `collector_args` | collector | all |

**`data`**: list of entries; an entry is one or more filterbank paths (several are merged into one
observation) optionally followed by the coherent DM (number), for example
`[["/data/obs1.fil", 56], ["/data/obs2_a.fil", "/data/obs2_b.fil", 0]]`. The coherent DM becomes every
pulsar's `cDM` (smearing) and is passed as `--cdm` to the PEASOUP search and the PulsarX parfolds and
candidate folds; PRESTO and DSPSR do not use it.

**`injection_args`**: `id` (tag added to the file name), `ephem`, `python` (executable in the image),
`cache`, `offline`, `generate_fb` (`[mean, std]` or `[]`), `filterbank` (`random` or `injection_number` to
cycle through `data`), `gulp_size_GB`, `stats_samples`, `save_merged_fb`, `save_inj_fb`, `save_report`,
`save_pars` (later steps need the injected filterbank and the report).

**Parfold sections**: `mode` (for example `psrfold_fil`), `use_predictor` (fold with the injector's
predictor: `t2pred` for PulsarX, `polycos` for PRESTO; use for extreme orbits), `channel_mask`
(`"f1:f2,f3:f4"` MHz), `template`, `blocksize_plan` and `nbin_plan` (values per period range,
`period_cutoff_s`), `nsubint`, `mask` (PRESTO), `cmd` (extra arguments), `save_*`, `transfer_TMP`,
`delete_inj_fb`.

**`filtool_args`**: `cmd` (cleaning options), `filplan` (frequency downsampling, bit depth, output mean and
std), `tscrunch` (list of time-downsampling factors).

**`peasoup_args`**: `filtool` (search the filtool output or the raw injected data), `ddplan`
(`"tscrunch": [DM_low, DM_high, DM_step]`), `segment_plan` (`"n": {acc_start, acc_end, min_snr,
fftsize}`; `"1"` full observation, `"2"` halves, ...), `inj_DM` and `inj_DM_offsets` (search only the
injected DMs, offsets like `[-0.05, 0, 0.05]`; include 0), `birdie_list` (`"f:width,..."` Hz),
`channel_mask`, `cmd`, `save_*`, `delete_*`. Jobs per injection = DDplan entries × total segments, for
example 3 × (1 + 2) = 9.

**GPU memory.** `fftsize` is either a list of prime factors (the pipeline takes the smallest product of
them that covers the segment) or one fixed length. The pipeline does not check the GPU. PEASOUP keeps the
whole FFT on the GPU, about 45 bytes per point plus cuFFT working space, which grows for lengths with large
powers of 3 and 5. Long, finely sampled segments can therefore run out of GPU memory: a 2-h observation at
15 µs (4.7 × 10⁸ samples) failed on 23 GB cards with `[2, 3, 5]` (3⁵ × 5⁹) and with 2²⁹, and ran with the
fixed length 2¹⁷ × 3⁶ × 5 = 477,757,440 (21.4 GB). For segments above about 3 × 10⁸ samples, check the
cards' memory (`nvidia-smi`), test the length on one job first, and if it does not fit, use a fixed length
dominated by powers of 2, shorter segments, or time-downsampled data (`tscrunch`).

**PRESTO sections**: `presto_rfi_args` (rfifind time and options, zero-DM birdie detection);
`presto_search_args` (`rfi_cleaner`: `rfifind`, `filtool` or `none`; `mask`, `birdies`, `bary`, `ddplan`,
`segment_plan` with `zmax`, `wmax`, `sigma`, `numharm`; `inj_DM`, `inj_DM_offsets`; extra options for
`prepdata`, `realfft`, `zapbirds`, `accelsearch`); `presto_sift_args` (`sigma_threshold`,
`c_pow_threshold`, `harm_pow_cutoff`, `r_err`, period limits, `min_num_DMs`, `low_DM_cutoff`, birdies;
lower `harm_pow_cutoff`, or set it to 0, for faint pulsars); `presto_candfold_args` (`fold_mode`:
`filterbank`, `filtool` or `dat`; `bary`, `mask`, extra `prepfold` options, products).

**`candidate_matcher_args`**: `mode` (`PEASOUP` or `PRESTO`), `harmonics` [`[0.5, 1, 2]`], `dm_level`
(DM window ends where the S/N falls to this fraction of the peak [0.3]), `max_folds` (matches folded per
pulsar, brightest first), `fold_all` (also fold the brightest unmatched candidates). A candidate matches if
its frequency or a harmonic lies within the pulsar's expected frequency range over the segment (orbit and
accelerations included) plus two Fourier bins, and its DM is inside the window.

**`classifier_args`**: `tag`, `python`, `script`, `model_dir`, `file_type`.

**`collector_args`**: one boolean per step (`pulsarx_parfold`, `presto_parfold`, `dspsr_parfold`,
`peasoup_search`, `pulsarx_candfold`, `presto_search`, `presto_candfold`, `sifter`, `classifier`). Set
`true` exactly for the steps that run in `main.nf`.

### Pipeline outputs

```text
output_dir/
  inj_000001/
    <beam>_<id>_<injection_id>.fil      injected filterbank
    report_<injection_id>_<seed>.json   injection report
    inj_pulsars/                        ephemerides, predictors; PULSARX/, PRESTO/, DSPSR/ parfolds
    processing/                         FILTOOL/, PEASOUP/, PRESTO/ intermediate products
    inj_cands/                          PEASOUP/<segment>/, PRESTO/<segment>/ candidate folds, match summaries
  RESULTS_<run>.csv                     one row per injected pulsar
  FILTERBANKS_<run>.txt
```

- Matches: `processing/<MODE>/inj_XXXXXX_<MODE>_matches.csv` and the summary
  `inj_cands/<MODE>/inj_XXXXXX_<MODE>_matches.json` (harmonic, frequency offset in bins, DM offset, DM
  window).
- `RESULTS_<run>.csv` column prefixes: `INJ_` injected parameters; `PAR_`, `PRESTO_PAR_`, `DSPSR_PAR_`
  parfold results; `PEA_`, `PRESTO_` best matched search candidate; `CAND_`, `PRESTO_CAND_` its fold and
  classifier scores.
- Recovery fraction = injections with a match / all injections, binned in injected parameters.

### Checking candidates yourself

The matcher converts frames internally. For your own comparisons:

- Injected `PEPOCH`, `T0`, `ACCEPOCH` (frame `bary`) are barycentric TDB epochs; the data start in that
  frame differs from the header `tstart` (UTC at the telescope) by TT−UTC (about 69 s) plus the Rømer
  delay (up to about 500 s).
- PEASOUP (and PRESTO unless `bary: true`) search topocentric data: f_topo = f_bary (1 − v_r/c), with v_r
  the telescope's velocity away from the pulsar (up to 30 km/s, several Fourier bins for fast pulsars).
- `processing/PEASOUP/inj_XXXXXX_PEASOUP_candidates.csv` columns: `period` (s), `dm`, `acc` (m/s²),
  `snr`, `xml_id`, `pepoch` (segment reference epoch, topocentric MJD), `n_samples`, `T` (FFT length, s),
  `XML`, `tscrunch`, `downsample`, `segment` (`i_n`), `fold_id`, `F_match` (1/`period`). Frequencies refer
  to `pepoch`. The `_matches.csv` adds `PSR_ID`, `harmonic`, `dF0_bins`, `dDM`.

Both corrections, for a topocentric MJD `t_topo` and a barycentric frequency `f_bary` (run from the
repository root):

```python
import json
from injector.setup_manager import SetupManager
report, fil = "report_<injection_id>_<seed>.json", "<injected>.fil"
ephem = json.load(open(report))["injection_report"]["ephem"]
models = SetupManager(report, fil, ephem, generate=True, stats_samples=0).pulsar_models
obs = models[0].obs
t_bary = obs.topo2bary([t_topo])[0]
v_r = obs.earth_radial_velocity([t_topo])[0]
f_topo = f_bary * (1 - v_r / 299792458.0)
```

### Turning steps on and off

Comment out the step's line in `main.nf` with `//`, then **pass the removed step's input channel to
whatever read its output**. Then set that step to `false` in `collector_args`.

- Only PulsarX parfold in INJECTOR: comment out `presto_parfold` and `dspsr_parfold`, and change
  `emit: done = inj_px.mix(inj_presto, inj_dspsr)` to `emit: done = inj_px`.
- No classifier: comment it out and use `inj_output = collapse_tag(inj_fold_cand)`.
- No candidate folding: comment out candfold and classifier, use
  `inj_output = collapse_tag(inj_tag_match)`.
- No filtool (PEASOUP): comment it out, use `inj_peasoup = PEASOUP(inj_pulsars)` and set
  `"filtool": false` in `peasoup_args`.
- No parfold: just comment it out (nothing reads it), but see the deletion issue in *Known issues*.
- A channel from `expand_plan` carries `(number, tag)`; it must pass `collapse_tag` before a step that
  expects a plain number (for example the collector).
- Unused `include` lines can stay. `-resume` reuses finished jobs of unchanged steps.

### Disk space

- Injected filterbanks are full size: 250 beams of 7.6 GB are 1.9 TB. Plan `delete_inj_fb` and
  `batch_size` with the user before launching, and estimate the peak disk use.
- `delete_inj_fb` and `delete_filtool_fb` in the search and
  candidate-folding sections delete them after the last job that needs them (tracked by marker files).
- `batch_size` limits how many exist at once.
- `save_inj_fb: false` keeps them out of `output_dir` (only if no later step needs them).
- Delete `work/` after a successful run (with approval).

### Cost of a run

Measured per injection in one campaign (10 pulsars, 12 PEASOUP searches per injection): injection 6.2 CPU-h,
PulsarX parfolds 5.4 CPU-h, PulsarX candidate folds 10.6 CPU-h, PEASOUP 4.9 GPU-h. An estimate from the
injection cost alone was 4 times too low. Parfolds scale with pulsars, candidate folds with pulsars ×
segments × `max_folds` (more with `fold_all`), searches with DDplan entries × segments. Time the
one-injection test and scale from it.

## Multi-epoch injections

- Fixed pulsar values, or a fixed `seed` per pulsar entry, so every injection gets the same pulsars.
- `PEPOCH`, `T0`, `POSEPOCH` as MJDs, not fractions, so spin and orbit continue between epochs; the same
  Solar System ephemeris for all epochs.
- PINT pulsars: the par file needs `TZRMJD` (normally with `TZRSITE`, `TZRFRQ`), or PINT picks a new phase
  reference per observation and the epochs are not phase-connected (INCEPTION prints a warning). For a
  made-up pulsar use e.g. `TZRMJD` = `PEPOCH`, `TZRSITE @`, `TZRFRQ 0`. Extra plan terms follow the point above.
- One `data` entry per epoch, `filterbank: "injection_number"`, `start` 1, `end` = number of epochs.
  Never reuse injection numbers within a campaign (the number sets the noise and the draws).
- Standalone injector with `--generate_fb`: a different `global_seed` per epoch, or every epoch gets the
  same noise.
- Synthetic epochs need a header with the right `tstart`; with `--generate_fb` only the input's header and
  length are used, so write a header-only (sparse) copy with a new `tstart`:

```python
import os, struct
src, dst, new_tstart = "template.fil", "epoch_02.fil", 60530.25
raw = open(src, "rb").read(65536)
end = raw.index(b"HEADER_END") + len(b"HEADER_END")
i = raw.index(b"tstart") + len(b"tstart")
header = raw[:i] + struct.pack("<d", new_tstart) + raw[i + 8:end]
with open(dst, "wb") as f:
    f.write(header)
    f.truncate(len(header) + os.path.getsize(src) - end)
```

- Check phase connection by folding every epoch with the same par file.

## Timing injected data

PINT injections are timing-grade from version 1.0.14. Measured: three pulsars at DM 200 in the ecliptic
(isolated; DM1+FD+NE_SW; ELL1), seven epochs over a year, MeerKAT L and UHF, DSPSR folds, per-channel PSRCHIVE
TOAs, PINT and TEMPO2 residuals against the unchanged par files without any fit: white at the 0.7–0.9 µs TOA
precision, reduced χ² 0.95–1.03, no annual, frequency or orbital structure, folded pulse at the model phase
within 0.2 µs. Directly against PINT, without folding: about 0.2 ns rms, at most 1 ns in every channel, from
110–190 MHz to 704–4032 MHz, 64–4096 channels; made-up DD binary (4.3-h orbit, e = 0.6, OMDOT, GAMMA, PBDOT,
Shapiro) 0.5 ns rms, 2.7 ns max; the same binary at 12 epochs over a year, timed end to end: reduced χ² 1.02
over 46,080 TOAs, no annual term, absolute offset −0.004 ± 0.026 µs, epoch means 0.13 µs rms against 0.09 µs
errors; Shapiro delay recovered with scale 1.003 ± 0.002. `EPHEM`, `CLK`, `PLANET_SHAPIRO` followed, also offline. Samples are evaluated
at their centre (`tstart` = start of the first sample).

Injection side, for timing-grade work:

- Inject from the par file the user will time with (`TZRMJD`, `TZRSITE`, `TZRFRQ` for several epochs).
  INCEPTION follows its `EPHEM` and `CLK`; recommend `EPHEM DE440` for new par files (DE421 against DE440:
  Earth about 116 km apart relative to the barycentre, constant offset up to 0.4 ms; built-in ephemeris
  about 10 µs).
- `--ephem` with the kernel named in `EPHEM`; PINT clock files current (one online run with `--cache`). A
  `TT(BIPMxxxx)` whose clock file PINT cannot load falls back to PINT's default realisation (about 1 ns
  apart) with a WARNING; to match exactly, put `tai2tt_bipmxxxx.clk` into the cache. `out/<ID>.par` keeps
  the `CLK` line (the user's timing software may have that file).
- A par file, not a ready-made `.polycos` file (used as is, at its own frequency only, not checked).
- Check the log's accuracy line (largest deviation from PINT across the band).
- Smearing within channels only with `DM_smear` `approx`/`exact` (default `off`), from |DM − cDM|; `cDM` has
  no effect otherwise. Coherently dedispersed data: `cDM` = DM; incoherent filterbanks: `DM_smear: "exact"`.
- Effects that move TOAs as for a real pulsar (useful to study their effect on timing): scattering,
  frequency-dependent (2D) profiles against one standard profile, `micro_structure` (experimental).
- Standalone injector with `--generate_fb`: a different `global_seed` per epoch.
- Fold and time with `out/<ID>.par`.

Fold and timing side, each worth 0.3–34 µs:

- Fold with a TEMPO2 predictor with enough coefficients:
  `tempo2 -f out/<ID>.par -pred "<site> <MJD1> <MJD2> <f_low> <f_high> 12 8 1800"`, DSPSR `-P`. DSPSR's
  own (2 frequency coefficients, 1-h segments) is off by µs. Compact or eccentric orbits: check the
  predictor against the par file, shorten the segments (last number, s) if needed (4.3-h orbit, e = 0.6:
  1800 s off by up to 10 µs, 300 s a few ns). PSRCHIVE config:
  `Predictor::default = tempo2`,
  `Predictor::policy = default`.
- Binaries in DSPSR: small blocks (`-U 1`); the default block (about a minute here) advances the phase
  linearly with the period at the block start (34 µs for a 0.3-d orbit).
- Binaries in PSRCHIVE: do not add sub-integrations (`pam -T`) before the channels are dedispersed; make
  TOAs from short sub-integrations. DSPSR folds all channels with one phase, so outer channels of a fast
  binary drift within a sub-integration (up to 1.3 µs TOA bias in the lowest L channels at periastron,
  2-min sub-ints); shorter sub-integrations reduce it.
- Standard profile of a binary: `psradd -ip` (add without re-alignment); the default re-aligns with
  psradd's own phase model, off by tens of bins for a fast orbit.
- Sub-band TOAs: PSRCHIVE dedisperses inside a sub-band without the Earth's Doppler shift (annual bias
  about 0.3 µs lowest L, 0.8 µs lowest UHF sub-band at DM 200, 8 channels each). Use per-channel TOAs or
  `Dispersion::barycentric_correction = 1`; also `WeightedFrequency::round_to_kHz = 0` (kHz rounding, 1.3 µs).
- TEMPO2 applies a tropospheric delay and `NE_SW` 4 when the par file has none; PINT and the injection
  apply neither unless asked. Add `CORRECT_TROPOSPHERE N` and an explicit `NE_SW` to par files used with
  TEMPO2.

## Follow-up mode (task E)

A trick that turns the validator pipelines into a tool for real data: the plan lists pulsars that are in
the data, or might be, each with `"SNR": 0`. Nothing is injected (data pass through unchanged when every
pulsar has `SNR: 0`), but ephemerides are written, so the pipeline folds the data at them (parfold),
searches as usual, matches search candidates to them and folds the matches.

- Known pulsars: `pint_polycos` par file (any timing or binary model).
- Unconfirmed candidates: in-built model, e.g. `P0`, `DM`, `A0` = 100 m/s², with `PEPOCH` at the
  candidate's reference epoch and `ACCEPOCH: 0`; a constant acceleration only holds within one
  observation.
- New pulsars: `candidate_matcher_args.fold_all: true` also folds the brightest unmatched candidates of
  every segment, so the run doubles as a normal search.

Ask: which pulsars or candidates (par files, or period, DM, acceleration and epoch), which observations,
which pipeline, and whether unmatched candidates should be folded too.

Settings: `data` = the observations, `filterbank: "injection_number"`, `start` 1, `end` = number of data
entries (each processed once); `generate_fb: []` (otherwise the data are replaced by noise); `inj_DM: true`
for a fast search around the listed DMs, but not for `pint_polycos` pulsars, whose plan `DM` is an offset
(see *Known issues*); each observation is copied to the output folder, so set `delete_inj_fb`.
Example: `injector/example_plans/09_followup.json`.

## Nullsar (task F)

Removes bright known pulsars by injecting a calibrated negative copy. Typical removal 90–99% of the
signal; for S/N above about 500, run it again on its output.

**Stages per tag:** `setup` → `filtool` → fold INIT (PulsarX at the ephemeris, full channel resolution) →
zap INIT (fit model: Gaussian components, per-channel amplitude, residual spin terms, light curve; inject
the negative model half a period away = calibration) → fold OPTIMISE → zap NULL (apply fitted amplitude,
phase, DM correction and per-subint light-curve correction; inject at the pulsar's phase) → fold CONFIRM.

**Inputs:**

- `nullsar/pipeline_general/configs/filterbanks.txt`: one line per file, `<tag> <path>`; same tag = one
  observation (merged).
- Par files (or PulsarX `.candfile`s) listed in `par_files`. They must be accurate at the epoch (frequency
  searches are off when folding).

**`nextflow.config`** (`nullsar/pipeline_general/`): `output_dir`, `pipeline_code` (path to
`nullsar/code`), `tmp_dir`, `filterbanks`, `pipeline_config` (`configs/null_parms.json`), `nullsar_image`
(Python), `pulsarx_image`. Cluster configs for several machines are in the folder (`hercules.config`,
`oz.config`, ...); copy one and adapt it as in *Pipeline set-up* (labels `setup`, `filtool`, `nullsar_fold`,
`nullsar_zap`).

**`configs/null_parms.json`**: `filtool_args.cmd`; in `nullsar`: `filtool` (run filtool first),
`SNR_limit`, `calibration_offset` [0.5], `dm_fit_sigma`, `dm_fit_min_chans`, `profile_components`
(`{"default": 2, "<psr>": "auto"}`), `par_files`, `injection` (`ephem`, `python`, `cache`, `offline`,
`gulp_size_GB`, `stats_samples`), `fold_pars` (`mode`, `channel_mask`, `template`, `blocksize`
(`"auto"` uses at least `min_block_pulses` pulses per block), `n_subint`, `cmd`), `delete_filtool`,
`delete_processing`.

**Run:**

```bash
cd nullsar/pipeline_general
nextflow run main.nf -c <cluster>.config
```

**Outputs** in `output_dir/<tag>/`: `01_FILES/` (input list, filtool output), `02_INIT/` (`FOLDS/`,
`MODELS/`, `INIT_fold_params.json`), `03_OPT/` (`FOLDS/`, `MODELS/`, **`*_NULL.fil`** = data with the
pulsars removed, ready to search), `04_CONFIRM/` (`FOLDS/`, `CONFIRM_<psr>.png`).

**Checks to report to the user:**

- `02_INIT/MODELS/model_<psr>.png`: residual panels should look like noise. Systematic profile residuals →
  more `profile_components` or `"auto"`.
- `03_OPT/MODELS/OPT_<psr>.png`: the inverted copy should sit about 0.5 turns from the pulse with similar
  shape and amplitude. A copy far from 0.5, or one that drifts or jumps in phase, means the model or the
  ephemeris does not match the data; report it to the user.
- CONFIRM fold: no pulse should remain. Leftovers: narrow +/− pattern (profile not captured), variation
  within sub-integrations (fast scintillation, time–frequency structure).
- Not removed by design: pulse-to-pulse variability. Very wide duty cycles make pulse and copy overlap.

## Verifying an injection

1. Fold the injected file with the written ephemeris, for example
   `psrfold_fil --template <template> --parfile out/<ID>.par -f out/<file>.fil` or
   `prepfold -noxwin -par out/<ID>.par out/<file>.fil`. For extreme binaries fold with the
   predictor (`t2pred` for PulsarX).
2. The pulse should be at one phase across time and frequency, with the expected profile.
3. Compare the folded S/N with the injected S/N for the chosen `SNR_def`. For isolated pulsars in clean
   Gaussian noise, PulsarX S/N / injected S/N was 0.94–1.18 (median 1.06) for 10 pulsars at S/N 22–146.
   Treat 0.9–1.2 as normal. Outside that range, check `SNR_def`, `DM_smear`, the channel mask and the noise
   statistics. Smearing, scattering, masks and RFI lower the ratio, as for a real pulsar.
4. Check the report: every pulsar, every resolved random value, the data statistics (`fb_mean`,
   `fb_sigma`).
5. Check the files: the injected file has the size of the input; the number of folds (`.ar`, PNG) equals
   the number of pulsars. A missing PNG next to an existing `.ar` points to a PulsarX issue (see *Known
   issues*).
6. For pure-noise campaigns, check that the noise is independent between injections: on about 20k samples
   × all channels of two files, the fraction of equal samples should be the chance value for rounded
   Gaussian noise, about 1/(2√π σ) (4% for σ = 7), and the correlation about 0.
7. Tell the user what you checked and show the fold image if one was made.

## Plan recipes

Isolated pulsar, realistic settings:

```json
{"psr_global": {"injection_id": "test", "global_seed": 1, "DM_smear": "exact", "DM_const": "historical"},
 "pulsars": [{"ID": "psr_A", "P0": 0.005, "DM": 50, "SNR": 25, "duty_cycle": 0.05, "spectrum": -1.6}]}
```

Binary with a predictor:

```json
{"ID": "bin_A", "P0": 0.0036, "DM": 80, "SNR": 25, "binary_period": 5.0, "M2": 0.2, "create_predictor": "t2pred"}
```

Cluster pulsar (own DM draw around the cluster DM, acceleration):

```json
{"ID": "msp_1", "P0": 0.00499, "DM": {"rng": "normal", "mean": 52.1, "sigma": 0.2}, "SNR": 30, "A0": 1.2}
```

Pure noise (task B): one pulsar with `"SNR": 0` plus `--generate_fb=<mean>,<std>`, with a new seed for
every file.

Follow-up (task E): every pulsar with `"SNR": 0`; known pulsars by `pint_polycos` par file, candidates by `P0`, `DM`, `A0`, `PEPOCH` and `ACCEPOCH: 0` (see *Follow-up mode*).

## Pitfalls

- `binary_period` is in **hours**; `scattering_time` in **ms**; `separation` in **arcmin**.
- `ACCEPOCH` is a fraction of the observation length counted from PEPOCH; `T0` is an MJD only.
- With `pint_polycos`, `DM` is an offset, not the DM.
- Draws in `psr_global` are shared by all pulsars of an injection; put per-pulsar draws in the pulsar
  entries.
- `DM_smear` defaults to `off`; set `exact` for realistic injections. `cDM` only enters the smearing, so
  with `DM_smear: "off"` it changes nothing.
- `DM_const` defaults to `exact`; set `historical` when the data are processed with TEMPO, PRESTO or
  PulsarX.
- `SNR` is at emission, before propagation; recovered S/N is lower with smearing and scattering.
- At most about 20 pulsars per filterbank, fewer for similar DMs, very high S/N or long periods
  (see *Pulsars per filterbank*).
- Inject into raw data, then clean (filtool, rfifind). Data that were cleaned with zapped samples replaced
  by constants get signal added to those constants, and block-normalising software (PulsarX) can turn them
  into a very bright channel.
- `--output` must exist. Relative paths in plans are resolved from the current directory; in the pipelines
  use absolute paths.
- In the pipelines, `collector_args` must match the steps in `main.nf`, and every path must be under a
  container bind (`-B`).
- Injected S/N is for the whole observation, folded. A search segment sees about
  S/N × √(T_seg/T_obs); Fourier searches see roughly the pulsed S/N (√(1 − DC) of the total for a top-hat
  pulse) minus harmonic-summing losses that depend on duty cycle and set-up. For a target search S/N, run
  a small pilot (PEASOUP with `inj_DM`) and scale from it.
- In the pipelines, the plan's `global_seed` is replaced per injection; a fixed `seed` inside a pulsar entry
  gives that pulsar the same draws in every injection. Never reuse injection numbers within a campaign.
- Compact binaries (orbital periods of a few hours or less): a TEMPO2 predictor that PulsarX builds from a
  par file is fitted beyond the observation, where the orbit or Taylor expansion can diverge. Fold with
  INCEPTION's own predictor (`create_predictor: "t2pred"`, `use_predictor: true`) and check that the
  predictor file exists.
- Keep pulsar IDs short (12 characters or less) for PulsarX folds (see *Known issues*).
- Every par file used for injection or folding needs `UNITS TDB` (see *Known issues*).

## Known issues

As of INCEPTION version 1.0.15.

| Issue | Effect | Work-around |
|---|---|---|
| In-built model not timing-grade (µs or better across epochs) | No Solar System Shapiro delay or topocentric TT–TDB term (up to tens of µs), TT(TAI) without observatory clock corrections while its par files state `CLK TT(BIPM)` (about 27 µs), fixed `EPHEM DE421` in par files (old TEMPO crashes with newer ones; against `--ephem de440.bsp` a constant fold offset up to 0.4 ms, changing < 1 µs per year), no solar-wind dispersion, DM delay without the Earth's Doppler shift (about 1e-4 of the delay) | Negligible for searches and folding. For timing, inject from a par file (`pint_polycos`), which is timing-grade. |
| `.polycos` files written by older INCEPTION versions | PINT wrote wrong reference phases for times before the par file's `TZRMJD`, so injections with such a file jump in phase every few minutes | Regenerate the polycos from the par file. |
| `inj_DM` with PINT pulsars | The search uses the plan's `DM`, which is an offset for `pint_polycos` pulsars | Use the full DM plan, or give the DMs explicitly. |
| Candidate folding's `delete_inj_fb` waits for the parfold | If the parfold is removed from `main.nf`, injected filterbanks are never deleted | Delete them yourself (with approval) or keep the parfold. |
| The injector does not create `--output` | It fails if the directory is missing | `mkdir -p` first. |
| The candidate matcher does not check acceleration | For compact binaries the expected frequency range spans many bins, so at low thresholds (about 10,000 candidates per segment at S/N 6) most matches are noise | Judge compact-binary matches by their folds and S/N; tell the user the match flag alone is unreliable there. |
| PulsarX writes no PNG, without an error, when the PNG path is longer than 90 characters | Parfolds give `.ar` and `.cands` but no PNG. The path is `<work_dir>/_<ID>/<ID>_<MJD>_<beam>_00001.png`, about len(work_dir) + 2 len(ID) + 40 characters, which passes 90 in a Nextflow scratch directory for IDs longer than about 12 characters | Keep pulsar IDs to 12 characters or less, or use a short `tmp_dir`. Check that the number of PNGs equals the number of pulsars after the test run. |
| PulsarX writes fixed files `/tmp/stdout.txt` and `/tmp/stderr.txt` when it builds a predictor from a par file | Parallel folds sharing `/tmp` can overwrite each other's files | Give each fold a private `/tmp` (SLURM `job_container/tmpfs`, or bind a per-job directory to `/tmp`). |
| Rarely (about 1 in 200 folds) a parfold PNG is missing although the fold succeeded | Missing PNG, data products fine | Re-run that fold if the PNG is needed. |
| Par files without a `UNITS` line are read as TCB by TEMPO2 and PulsarX | PEPOCH shifts by about 23 s and F0 by F1 × 23 s; bad folds for large F1 (for example compact orbits written as spin derivatives) | Add `UNITS TDB` to every par file used for folding or with `pint_polycos`. INCEPTION's own par files contain it; PINT injections copy the input par file as it is. |

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| One channel much brighter in the fold | Input was cleaned before injection (see *Pitfalls*). Inject before cleaning or mask the channel. |
| Compact binary folds with wiggles in PulsarX | PulsarX's internal predictor cannot follow the orbit. `create_predictor: "t2pred"` and `use_predictor: true`. |
| Recovered DM off by about 0.01% | `DM_const`; use `historical`. |
| Job hangs after printing the filterbank statistics | Waiting for downloads. Use `--cache` and `--offline`. |
| `--offline` fails at once for PINT | Clock files not cached. Run once online with `--cache`. |
| Phase jumps every few minutes in a PINT injection | A `.polycos` file from an older INCEPTION version (see *Known issues*); regenerate it from the par file. |
| Annual or frequency-dependent trend in timing residuals of a PINT injection | `.polycos` covers only the band centre (user-supplied, INCEPTION 1.0.11 or earlier, or `pint_nfreq: 1`); regenerate from the par file with `pint_nfreq` above 1. |
| Constant offset of about 27.7 µs, or of half a sampling time, in timing residuals of a PINT injection | INCEPTION 1.0.12 or earlier: polycos ignored the par file's `CLK` (27.67 µs for `TT(TAI)`) and samples were evaluated at their start. Fixed in 1.0.13; reinject. |
| PINT injections at different epochs not phase-connected | Par file has no `TZRMJD` (warning in the log); add `TZRMJD`, `TZRSITE`, `TZRFRQ` (see *Multi-epoch injections*). |
| `create_parfile` rejected | Must be `"par"` or `"pulsarx"`. |
| Recovered S/N differs from injected | `SNR_def` against the folding software; smearing, scattering, masks and RFI lower it. |
| Nextflow process fails | Read `.command.err` and `.command.log` in the work directory named in the error. Path missing in a container: add a `-B` bind. |
| PEASOUP job fails with "No xml file produced" | PEASOUP did not finish; look in `.command.log` for the cause. `CUDA failed with error: out of memory`: FFT length too large for the GPU (see *GPU memory* under `pipeline_configs.json`). |
| Matcher stops with "Ephemeris file not found" | The report points to a `.polycos` file that no longer exists, typically a report from an older INCEPTION version, which pointed to the injection job's scratch copy. Point the report's `pint_polycos` to the copy in `inj_pulsars/`, or rerun the injection. |
| Collector fails or columns empty | `collector_args` has a step that is not in `main.nf`. |
| Injected filterbanks not deleted | Parfold removed (see *Known issues*). |

## Contact and citing

For problems this guide does not solve, suggest the user contacts Robert Senzel
(rsenzel@mpifr-bonn.mpg.de), the INCEPTION developer, with the plan, the command, the log and the
injection report.

Ask the user to cite the INCEPTION paper (Senzel et al., in prep.) and the software the pipelines call
(PulsarX, PRESTO, DSPSR/PSRCHIVE, PEASOUP, PINT, Astropy, Nextflow).
