import sys
import json
import subprocess
from collections import namedtuple



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
    

def create_DDplan(ddplan):
    DMRange = namedtuple('DMRange', ['low_dm', 'high_dm', 'dm_step', 'tscrunch'])

    segments = []
    for line in ddplan.splitlines():
        low_dm, high_dm, dm_step, tscrunch = list(map(float, line.split()[:4]))
        segments.append(DMRange(low_dm, high_dm, dm_step, tscrunch))

    return list(sorted(segments, key=lambda x: x.tscrunch))


