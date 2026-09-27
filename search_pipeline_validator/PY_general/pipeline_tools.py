import os
import sys
import glob
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

def glob_psr(directory, psr_id, suffix):
    psr_id = glob.escape(str(psr_id))
    return sorted(glob.glob(f'{directory}/{psr_id}{suffix}') +
                  glob.glob(f'{directory}/{psr_id}_*{suffix}'))

def string2seed(s):
    hash_object = hashlib.sha256(s.encode())
    hash_int = int(hash_object.hexdigest(), 16)
    return hash_int % (10**12)


def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + str(output))


def presto_rfi_cleaner(processing_args):
    s_args = processing_args['presto_search_args']
    cleaner = s_args.get('rfi_cleaner', 'rfifind')
    if cleaner not in ('rfifind', 'filtool'):
        sys.exit(f"presto_search_args.rfi_cleaner must be 'rfifind' or 'filtool', not '{cleaner}'.")

    if cleaner == 'filtool':
        conflicts = []
        if s_args.get('mask'):
            conflicts.append('presto_search_args.mask')
        if processing_args.get('presto_candfold_args', {}).get('mask'):
            conflicts.append('presto_candfold_args.mask')
        if s_args.get('birdies') == 'rfifind':
            conflicts.append('presto_search_args.birdies')
        if processing_args.get('presto_parfold_args', {}).get('mask') == 'rfifind':
            conflicts.append('presto_parfold_args.mask')
        if conflicts:
            sys.exit(f"rfi_cleaner is 'filtool', so {', '.join(conflicts)} must not use rfifind or a mask "
                     f"- only one RFI cleaner can be used.")
    return cleaner

def filtool_filterbank(results_dir, inj_id, processing_args):
    tscrunch = processing_args['filtool_args']['tscrunch']
    if 1 not in tscrunch:
        sys.exit('filtool_args.tscrunch must include 1: PRESTO dedisperses the full-resolution filtool output.')

    data = glob.glob(f"{results_dir}/processing/FILTOOL/*_{inj_id}_FILTOOL_0{tscrunch.index(1) + 1}.fil")
    if not data:
        sys.exit('No filtool-cleaned filterbank found - filtool_args.save_filtool_fb must be true.')
    return data[0]

def parse_process_tag(process_tag):
    splits = process_tag.split('_')
    downsample = int(splits[3])

    if 'SEG' in splits:
        return downsample, int(splits[5]), int(splits[6])
    return downsample, None, None


def build_dm_list(ddplan, downsample, inj_DM, injection_report):
    if inj_DM:
        return sorted(psr['DM'] for psr in injection_report['pulsars'])

    low, high, step = ddplan[str(downsample)]
    n_trial = int(round((high - low) / step))
    return list(np.linspace(low, high, n_trial, endpoint=False))


def segment_samples(n_total, seg_i, seg_n):
    start_sample = int(np.floor(seg_i * n_total / seg_n))
    end_sample = int(np.floor((seg_i + 1) * n_total / seg_n))
    return start_sample, end_sample - start_sample


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