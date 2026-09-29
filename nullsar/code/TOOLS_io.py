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
            if (len(columns) < 2) or columns[0].startswith('#') or (columns[0] == 'C'):
                continue
            data.setdefault(columns[0], columns[1])
    return data


def par_float(value):
    return float(str(value).replace('D', 'E').replace('d', 'e'))


def par_period(params):
    if 'F0' in params:
        return 1/par_float(params['F0'])
    return par_float(params['P0'])


def set_par_value(parfile, new_parfile, key, value):
    with open(parfile, "r") as f:
        lines = f.readlines()

    replaced = False
    for i, line in enumerate(lines):
        columns = line.split()
        if columns and (columns[0] == key):
            lines[i] = ' '.join([key, str(value), *columns[2:]]) + '\n'
            replaced = True

    if not replaced:
        if lines and not lines[-1].endswith('\n'):
            lines[-1] += '\n'
        lines.append(f'{key} {value}\n')

    with open(new_parfile, "w") as f:
        f.writelines(lines)


def add_cmd_args(cmd, args):
    extra = args.get('cmd', '')
    if not isinstance(extra, str) or 'cmd_flags' in args:
        sys.exit(f"'cmd' must be one string of arguments, e.g. \"--nbin 64 -v\", and 'cmd_flags' is no longer used: {args}")
    return f'{cmd} {extra}'


def parse_cand_file(candifle):
    with open(candifle, 'r') as f:
        data = [line.strip().split() for line in f if line.strip()]
    return dict(zip(data[0], data[1]))