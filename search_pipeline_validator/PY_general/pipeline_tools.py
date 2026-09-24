import os
import sys
import json
import hashlib
import subprocess
import numpy as np



def parse_JSON(json_file):
    try:
        with open(json_file, 'r') as file:
            pars = json.load(file)
    except FileNotFoundError:
        sys.exit(f'Unable to find {json_file}.')
    except json.JSONDecodeError:
        sys.exit(f'Unable to parse {json_file} using JSON.')
    else:
        return pars
    
def rsync(source, destination, shell=True):
    try:
        subprocess.run(f'rsync -Pav {source} {destination}', shell=shell, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        subprocess.run(f'cp -av {source} {destination}', shell=shell)
    
def next_fast_len(n, primes=[2, 3, 5]):
    """Return the smallest integer >= n with prime factors only in `primes`."""

    def is_smooth(x):
        for p in primes:
            while x % p == 0:
                x //= p
        return x == 1

    m = n
    while not is_smooth(m):
        m += 1

    return m

def string2seed(s):
    hash_object = hashlib.sha256(s.encode())
    hash_int = int(hash_object.hexdigest(), 16)
    return hash_int % (10**12)


def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + str(output))


def build_dm_trials(ddplan, inj_DM, injection_report):
    """Build the (DM, downsample) dedispersion job list. `ddplan` is keyed by
    downsample factor (as a string, e.g. peasoup's ddplan convention), each value
    a [low, high, step] DM sweep for that downsample - so each downsample gets
    its own DM grid. If inj_DM is True, every downsample uses the exact injected
    DM values instead of sweeping its [low, high, step] range."""
    trials = []
    for ds_str, dm_range in ddplan.items():
        ds = int(ds_str)
        if inj_DM:
            DM_values = sorted(psr['DM'] for psr in injection_report['pulsars'])
        else:
            low, high, step = dm_range
            n_trial = int(round((high - low) / step))
            DM_values = np.linspace(low, high, n_trial, endpoint=False)
        trials.extend((float(dm), ds) for dm in DM_values)
    return trials


def batch_trials(trials, batch_size):
    """Split a flat trial list into fixed-size batches, deterministically (same
    trials + batch_size always yields the same batches, so every pipeline stage
    can independently recompute which trials a given batch tag covers)."""
    return [trials[i:i + batch_size] for i in range(0, len(trials), batch_size)]