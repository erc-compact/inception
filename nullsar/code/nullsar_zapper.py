import json
import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path

from ar_processor import ARProcessor
from nullsar_tools import parse_cand_file, parse_par_file, parse_JSON, rsync, fit_time_phase, fit_phase_offset, scale_freq_phase

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from injector import SCRIPT_inject_pulsars 
from injector.io_tools import merge_filterbanks, print_exe, FilterbankReader


class InjectorProcess:
    def __init__(self, tag, processing_args, out_dir, work_dir, mode='INIT'):

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.processing_args_path = processing_args
        self.processing_args = parse_JSON(processing_args)['nullsar']

        self.tag = tag
        self.mode = mode
        self.processing_dir = f'{self.out_dir}/PROCESSING/{self.tag}'


    def get_data(self):
        files_dir = f'{self.processing_dir}/01_FILES'

        if self.processing_args.get('filtool', False):
            self.data = glob.glob(f'{files_dir}/*FILTOOL*.fil')[0]
        else:
            with open(f'{files_dir}/files.txt') as f:
                self.data = [line.strip() for line in f if line.strip()]

    def extract_archive(self):
        par_files = self.processing_args['par_files']
        files_dir = f'{self.processing_dir}/01_FILES/NULLSAR'
        ar_path = f"{files_dir}/INIT_fold_params.json"
        self.ar_data = {}

        for par_file in par_files:
            psr_ID = Path(par_file).stem
            if self.mode == 'INIT':
                fits_path = f'{files_dir}/{psr_ID}_mode_INIT.fits'
            if self.mode == 'NULL':
                fits_path = f'{files_dir}/{psr_ID}_mode_OPTIMISE.fits'

            archive = ARProcessor(fits_path, mode='load')
            self.parse_archive(psr_ID, archive, ar_path)

        if self.mode == "INIT":
            with open(ar_path, 'w') as file:
                json.dump(self.ar_data, file, indent=4)

    def parse_archive(self, psr_ID, archive, ar_path):
        files_dir = f"{self.processing_dir}/01_FILES/NULLSAR"

        fb = FilterbankReader(self.new_fb_path, load_fb_stats=(128, 6))
        obs_len = fb.dt * fb.n_samples

        profile_path = f"{files_dir}/profile_{psr_ID}.npy"        
        if self.mode == 'INIT':
            SNR = archive.get_SNR()
            freq_phase = archive.get_freq_phase()
            freq_phase_scaled = scale_freq_phase(freq_phase)

            time_phase = archive.get_time_phase()
            freq_deriv, phase_offset = fit_time_phase(time_phase, freq_phase, obs_len)

            np.save(profile_path, freq_phase_scaled)

        elif self.mode == "NULL":
            init_ar_data =  parse_JSON(ar_path)
            SNR = init_ar_data[psr_ID]['SNR']
            freq_deriv = init_ar_data[psr_ID]['FX']

            intensity_profile = archive.get_intensity_prof()
            freq_phase = np.load(profile_path)

            phase_offset, SNR_scale = fit_phase_offset(intensity_profile, freq_phase)
            SNR *= SNR_scale
            phase_offset += init_ar_data[psr_ID]['phase_offset']

        self.ar_data[psr_ID] = {"SNR": SNR,  
                                "phase_offset": phase_offset, 
                                "profile": profile_path,
                                "FX": freq_deriv}
            
    def create_injection_plan(self):
        injection_plan = {
            "psr_global": {
                "injection_id": self.mode,
                "global_seed": 404,
                "create_parfile": 0,
            },

            "pulsars": []
        }

        par_files = self.processing_args['par_files']
        for par_file in par_files:
            psr_ID = Path(par_file).stem
            params = parse_par_file(par_file)

            psr_dict = {
                "ID": psr_ID,
                "mode": "pint",
                "PEPOCH": 0.5,
                "phase_offset": self.ar_data[psr_ID]['phase_offset'],
                
                "P0_SNR": 1/float(params['F0']),
                "DM": float(params['DM']),
                "SNR": -self.ar_data[psr_ID]['SNR'],

                "profile": self.ar_data[psr_ID]['profile'],
                "polycos": par_file
            }

            for key, value in self.ar_data[psr_ID]['FX'].items():
                psr_dict[key] = value

            injection_plan['pulsars'].append(psr_dict)

        files_dir = f"{self.processing_dir}/01_FILES/NULLSAR"
        self.inject_file = f"{files_dir}/NULLSAR_inject_file_mode_{self.mode}.json"
        with open(self.inject_file, 'w') as file:
            json.dump(injection_plan, file, indent=4)

    def injector_setup(self):
        self.get_data()
        self.transfer_merge()
        self.extract_archive()
        self.create_injection_plan()

    def transfer_merge(self):
        if type(self.data) == list:
            fb_names = [Path(fb).stem for fb in self.data]

            prefix = os.path.commonprefix(fb_names)
            self.new_fb_path = f"{self.work_dir}/{prefix}_MERGED.fil"

            for data_product in self.data:
                rsync(data_product, self.work_dir)
            
            new_data_paths = [f'{self.work_dir}/{Path(fb).name}' for fb in self.data]
            new_data_paths.sort()
            merge_filterbanks(new_data_paths, self.new_fb_path)

        else:
            prefix = Path(self.data).stem 
            self.new_fb_path = f"{self.work_dir}/{prefix}.fil"

            rsync(self.data, self.new_fb_path)

        ephem = self.processing_args['injection']['ephem']
        if ephem != 'builtin':
            rsync(ephem, self.work_dir)
            self.ephem = f'./{Path(ephem).name}'

    def run_injector(self, ncpus):
        gulp_size = self.processing_args['injection']['gulp_size_GB']
        n_samples = self.processing_args['injection']['stats_samples']

        inputs = f"--signal={self.inject_file} --fb={self.new_fb_path} --ephem={self.ephem} --output={self.work_dir} --ncpu={ncpus} --gulp_size_GB={gulp_size} --stats_samples={n_samples}"
        cmd = f"{self.processing_args['injection']['python']} {SCRIPT_inject_pulsars.__file__} {inputs}"

        print_exe('starting injection...')
        subprocess.run(cmd, shell=True)
        print_exe('injection complete.')


    def transfer_products(self):
        injected_fb = f'{self.work_dir}/{Path(self.new_fb_path).stem}_{self.mode}.fil'
        if not glob.glob(injected_fb):
            sys.exit(1)

        results_dir = f'{self.processing_dir}/01_FILES/NULLSAR'
        os.makedirs(results_dir, exist_ok=True)

        rsync(injected_fb, results_dir)
        rsync(f'{self.work_dir}/*.polycos', results_dir)

        

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='COMPASS - preprocessing channel splitter',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='file tag')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--mode', metavar='int', type=str, required=False, default="init", help='folding mode (INIT, NULL)')
    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of cpus for injection')

    args = parser.parse_args()
    inj_exec = InjectorProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.mode)
    inj_exec.injector_setup()
    inj_exec.run_injector(args.ncpus)
    inj_exec.transfer_products()