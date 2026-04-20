import sys, os
import json
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

def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + str(output))
    

def rsync(source, destination, shell=True):
    try:
        subprocess.run(f'rsync -Pav {source} {destination}', shell=shell, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        subprocess.run(f'cp -av {source} {destination}', shell=shell)


def parse_par_file(parfile):
    data = {}
    with open(parfile, "r") as f:
        for line in f:
            columns = line.split()
            data[columns[0]] = columns[1]
    return data


def parse_cand_file(candifle):
    with open(candifle, 'r') as f:
        data = [line.strip().split() for line in f if line.strip()]
    return dict(zip(data[0], data[1]))