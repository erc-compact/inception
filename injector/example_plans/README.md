# Example injection plans

Run them from the repository root, for example

```bash
mkdir -p out
python3 injector/SCRIPT_inject_pulsars.py --signal=injector/example_plans/01_minimal.json --fb=obs.fil --output=out/ --ncpu=8
```

File paths inside the plans (`gain_map`, `profile`, `pint_polycos`) are relative to the directory the
injector is run from. Use absolute paths when running elsewhere. See the INCEPTION user guide for every
parameter. Inject about 10, and at most about 20, pulsars per filterbank; build larger populations
with the Nextflow pipelines, which run the same plan on many filterbanks.

| Plan | Shows |
|---|---|
| `01_minimal.json` | The smallest valid plan: an isolated pulsar with period, DM and S/N. |
| `02_binaries.json` | A relativistic double neutron star given by masses and inclination, a neutron star–white dwarf binary given by its companion mass, and an eclipsing spider given by its projected semi-major axis and an eclipse light curve (`eclipse.npy`). Writes `t2pred` predictors for folding. |
| `03_double_pulsar.json` | Pulsar B is made the companion of pulsar A with `double_pulsar`. |
| `04_cluster_accelerations.json` | Globular-cluster pulsars with DMs drawn around the cluster DM, their own spectral indices, shared smearing and dispersion constant in `psr_global`, line-of-sight acceleration and jerk, and PRESTO `z` units. |
| `05_custom_profiles.json` | Fixed Gaussian components, randomised components, and a profile from a file (`example_profile.npy`). |
| `06_propagation.json` | Scattering with a given timescale, scattering from the DM–scattering relation, and smearing in coherently dedispersed data (`cDM`). |
| `07_population.json` | A replicated population of ten pulsars with random draws (uniform, log-uniform, normal, weighted choice, orbital periods in units of the observation length). |
| `08_known_pulsar_pint.json` | A known pulsar injected from its timing ephemeris with PINT, and the same pulsar with a DM offset and an extra acceleration. Replace `path/to/known_psr.par`. Needs `pint-pulsar`. |
| `09_followup.json` | Follow-up mode: every pulsar has `SNR: 0`, so the data pass through unchanged while ephemerides are written for folding and candidate matching. Replace `path/to/candidate_1.par`. |
