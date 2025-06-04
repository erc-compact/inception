import json
import argparse
import subprocess
import numpy as np
from pathlib import Path

import pipeline_tools as inj_tools

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from injector import SCRIPT_inject_pulsars 
from injector.io_tools import merge_filterbanks, print_exe


class InjectorProcess:
    def __init__(self, processing_args, injection_plan, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.injection_plan_path = injection_plan
        self.injection_plan = inj_tools.parse_JSON(self.injection_plan_path)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number
        self.rng = np.random.default_rng(injection_number)
        
    def injector_setup(self):
        self.resolve_seed()
        self.choose_data()
        self.transfer_merge()

    def resolve_seed(self):
        self.seed = int(self.rng.integers(1e11, 1e12))

        self.seeded_inject_file = f'{self.work_dir}/{Path(self.injection_plan_path).name}'
        self.injection_plan['psr_global']['global_seed'] = self.seed
        with open(self.seeded_inject_file, 'w') as file:
            json.dump(self.injection_plan, file, indent=4)

    def choose_data(self):
        beams = self.processing_args['data']
        if self.processing_args['injection_args']['filterbank'] == 'random':
            filterbanks = beams[self.rng.integers(0, len(beams))]
        elif self.processing_args['injection_args']['filterbank'] == 'injection_number':
            filterbanks = beams[(self.injection_number-1) % len(beams)]

        self.data = filterbanks

    def transfer_merge(self):
        if type(self.data) == list:
            fb_names = [Path(fb).stem for fb in self.data]

            prefix = os.path.commonprefix(fb_names)
            self.new_fb_path = f"{self.work_dir}/{prefix}_MERGED_{self.processing_args['injection_args']['id']}.fil"

            for data_product in self.data:
                inj_tools.rsync(data_product, self.work_dir)
            
            new_data_paths = [f'{self.work_dir}/{Path(fb).name}' for fb in self.data]
            new_data_paths.sort()
            merge_filterbanks(new_data_paths, self.new_fb_path)

        else:
            prefix = Path(self.data).stem 
            self.new_fb_path = f"{self.work_dir}/{prefix}_{self.processing_args['injection_args']['id']}.fil"

            inj_tools.rsync(self.data, self.new_fb_path)

    def run_injector(self, ncpus):
        ephem = self.processing_args['injection_args']['ephem']
        if ephem != 'builtin':
            inj_tools.rsync(ephem, self.work_dir)
            ephem = f'./{Path(ephem).name}'
        gulp_size = self.processing_args['injection_args']['gulp_size_GB']
        n_samples = self.processing_args['injection_args']['stats_samples']

        inputs = f"--signal={self.seeded_inject_file} --fb={self.new_fb_path} --ephem={ephem} --output={self.work_dir} --ncpu={ncpus} --gulp_size_GB={gulp_size} --stats_samples={n_samples}"
        cmd = f"{self.processing_args['injection_args']['python']} {SCRIPT_inject_pulsars.__file__} {inputs}"

        print_exe('starting injection...')
        subprocess.run(cmd, shell=True)
        print_exe('injection complete.')


    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        par_dir = f'{results_dir}/inj_pulsars'

        inj_ID = self.injection_plan['psr_global']['injection_id']
        injected_fb = f'{self.work_dir}/{Path(self.new_fb_path).stem}_{inj_ID}.fil'
        injection_report = f'{self.work_dir}/report_{inj_ID}_{self.seed}.json'
        
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(par_dir, exist_ok=True)
        
        if self.processing_args['injection_args']['save_merged_fb']:
            inj_tools.rsync(self.new_fb_path, results_dir)
        if self.processing_args['injection_args']['save_inj_fb']:
            inj_tools.rsync(injected_fb, results_dir)
        if self.processing_args['injection_args']['save_report']:
            inj_tools.rsync(injection_report, results_dir)
        if self.processing_args['injection_args']['save_pars']:
            inj_tools.rsync(f'{self.work_dir}/*.par', par_dir)
        

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='INCEPTION - injector process',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_plan', metavar='file', required=True, help='JSON file with injection plan')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of cpus for injection')
    args = parser.parse_args()


    inj_exec = InjectorProcess(args.processing_args, args.injection_plan, args.out_dir, args.work_dir, args.injection_number)
    inj_exec.injector_setup()
    inj_exec.run_injector(args.ncpus)
    inj_exec.transfer_products()