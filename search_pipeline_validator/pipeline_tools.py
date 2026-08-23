import os
import sys
import json
import hashlib
import subprocess



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