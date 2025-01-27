import os
import json
import argparse
import subprocess
import numpy as np
from pathlib import Path

from pipeline_tools import PipelineTools
from inception.injector import SCRIPT_inject_pulsars 
from inception.injector.io_tools import merge_filterbanks, print_exe


class InjectorSetup(PipelineTools):
    def __init__(self, search_args, inject_file, data_dir, out_dir, injection_number):
        super().__init__(search_args)
        self.work_dir = os.getcwd()
        self.data_dir = data_dir
        self.out_dir = out_dir
        self.injection_number = injection_number
        self.rng = np.random.default_rng(injection_number)

        self.parse_injection_file(inject_file)
        self.resolve_seed()
        self.get_beam()
        self.rsync_merge_data_products()
    
    def parse_injection_file(self, inject_file):
        self.inject_file_path = inject_file
        self.injected_params = self.parse_JSON(inject_file)

    def resolve_seed(self):
        self.seed = int(self.rng.integers(1e11, 1e12))

        self.seeded_inject_file = f'{self.work_dir}/{Path(self.inject_file_path).name}'
        self.injected_params['psr_global']['global_seed'] = self.seed
        with open(self.seeded_inject_file, 'w') as file:
            json.dump(self.injected_params, file, indent=4)
        
    def get_beam(self): 
        beams = []
        for pointing in self.data['pointings']:
            for beam in pointing['beams']:
                beam_data = [pointing['id'], beam['name']]
                for filterbank in beam['data_products']:
                    beam_data.append(filterbank['filename'])
                beams.append(beam_data)

        self.beam_data = beams[self.rng.integers(0, len(beams))]
    
    def rsync_merge_data_products(self):
        pointing_id, inj_beam_name, *fb_names = self.beam_data
        self.merged_fb = f'{self.work_dir}/{pointing_id}_{inj_beam_name}_merged_{self.data_ID}.fil'

        data_paths = [f'{self.data_dir}/{pointing_id}/{inj_beam_name}/{fb_name}' for fb_name in fb_names]
        for data_product in data_paths:
            cmd = f"rsync -PavL {data_product} {self.work_dir}"
            subprocess.run(cmd, shell=True)
            print(cmd)

        new_data_paths = [f'{self.work_dir}/{fb_name}' for fb_name in fb_names]
        merge_filterbanks(new_data_paths, self.merged_fb)

    def run_injector(self, ephem, ncpus):
        subprocess.run(f"rsync -Pav {ephem} {self.work_dir}", shell=True)
        
        inputs = f"--signal={self.seeded_inject_file} --fb={self.merged_fb} --ephem=./de440.bsp --output={self.work_dir} --ncpu={ncpus}"
        cmd = f"python3 {SCRIPT_inject_pulsars.__file__} {inputs}"
        print_exe('starting injection')
        try:
            subprocess.run(cmd, shell=True)
        except:
            raise
        print_exe('injection complete')

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        par_dir = f'{results_dir}/inj_pulsars'
        inj_ID = self.injected_params['psr_global']['injection_id']
        injected_fb = f'{self.work_dir}/{Path(self.merged_fb).stem}_{inj_ID}.fil'
        injection_report = f'{self.work_dir}/report_{inj_ID}_{self.seed}.json'
        
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(par_dir, exist_ok=True)
        # subprocess.run(f"rsync -Pav {self.merged_fb} {results_dir}", shell=True)
        subprocess.run(f"rsync -PavL {injected_fb} {results_dir}", shell=True)
        subprocess.run(f"rsync -PavL {injection_report} {results_dir}", shell=True)
        subprocess.run(f"rsync -PavL {self.work_dir}/*.par {par_dir}", shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='offline injection pipeline setup',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_file', metavar='file', required=True, help='JSON file with injection plan')
    parser.add_argument('--data_dir', metavar='dir', required=True, help='directory of observation filterbanks')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--ephem', metavar='file', required=False, default='builtin', help='solar system ephemeris file')
    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of threads for injection')
    args = parser.parse_args()

    inj_setup = InjectorSetup(args.search_args, args.injection_file, args.data_dir, args.out_dir, args.injection_number)
    inj_setup.run_injector(args.ephem, args.ncpus)
    inj_setup.transfer_products()