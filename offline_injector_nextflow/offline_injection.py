import os
import sys
import json
import argparse
import subprocess
import numpy as np

from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parent.parent))
from inception.injector.io_tools import merge_filterbanks


class InjectorSetup:
    def __init__(self, search_args, inject_file, data_dir, work_dir):
        self.out = work_dir
        self.data_dir = data_dir

        args = self.parse_JSON(search_args)
        self.inject_file = self.parse_JSON(inject_file)
        self.processing_args, self.ID  = args['processing_args'], args['processing_id']
        self.beam_data = self.get_beam(args['data']['pointings'])
        
        self.injection_ID = self.inject_file['psr_global']['injection_id']
        self.data_ID = args['processing_id']

    def parse_JSON(self, json_file):
        try:
            with open(json_file, 'r') as file:
                pars = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {json_file}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {json_file} using JSON.')
        else:
            return pars
    
    def get_beam(self, pointings):
        beams = []
        for pointing in pointings:
            for beam in pointing['beams']:
                beam_data = [pointing['id'], beam['name']]
                for filterbank in beam['data_products']:
                    beam_data.append(filterbank['filename'])
                beams.append(beam_data)

        selected_beam = beams[np.random.randint(0,len(beams))]
        return selected_beam
    

    def rsync_merge_data_products(self):
        pointing_id, inj_beam_name, *fb_names = self.beam_data
        output_name = f'{self.out}/{pointing_id}_{inj_beam_name}_{self.data_ID}_merged.fil'

        data_paths = [f'{self.data_dir}/{pointing_id}/{inj_beam_name}/{fb_name}' for fb_name in fb_names]
        for data_product in data_paths:
            cmd = f"rsync -Pav {data_product} {self.out}"
            subprocess.run(cmd, shell=True)

        new_data_paths = [f'{self.out}/{fb_name}' for fb_name in fb_names]
        merge_filterbanks(new_data_paths, output_name)
        return output_name

    def run_injector(self, injection_file, ephem, ncpus):
        fb_path = self.rsync_merge_data_products()

        script_path = '/hercules/u/rsenzel/pulsar_inject/inception/injector'
        inputs = f"--signal={injection_file} --fb={fb_path} --ephem={ephem} --output={self.out} --ncpu={ncpus}"
        cmd = f"python3 {script_path}/SCRIPT_inject_pulsars.py {inputs}"

        subprocess.run(cmd, shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='offline injection pipeline setup',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_file', metavar='file', required=True, help='JSON file with injection plan')
    parser.add_argument('--data_dir', metavar='dir', required=True, help='directory of xml files')
    parser.add_argument('--work_dir', metavar='dir', required=True, help='work directory')
    parser.add_argument('--ephem', metavar='file', required=False, default='builtin', help='solar system ephemeris file')
    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of threads for injection')
    args = parser.parse_args()

    inj_setup = InjectorSetup(args.search_args, args.injection_file, args.data_dir, args.work_dir)
    inj_setup.run_injector(args.injection_file, args.ephem, args.ncpus)


