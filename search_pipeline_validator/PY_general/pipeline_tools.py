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


def parse_process_tag(process_tag):
    splits = process_tag.split('_')
    return int(splits[3]), int(splits[5]), int(splits[6])


def build_dm_list(ddplan, downsample, inj_DM, injection_report):
    if inj_DM:
        return sorted(psr['DM'] for psr in injection_report['pulsars'])

    low, high, step = ddplan[str(downsample)]
    n_trial = int(round((high - low) / step))
    return list(np.linspace(low, high, n_trial, endpoint=False))


def segment_samples(n_total, seg_i, seg_n, downsample):
    start_sample = int(np.floor(seg_i * n_total / seg_n))
    end_sample = int(np.floor((seg_i + 1) * n_total / seg_n))

    start_frac = start_sample / n_total
    numout = (end_sample - start_sample) // downsample
    return start_frac, numout


def add_cmd_args(cmd, args, skip_flags=(), skip_keys=()):
    for flag in args.get('cmd_flags', []):
        if flag in skip_flags:
            continue
        cmd += f' {flag}'

    for key, value in args.get('cmd', {}).items():
        if key in skip_keys:
            continue
        cmd += f' -{key} {value}'

    return cmd