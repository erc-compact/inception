import json
import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path

from TOOLS_ar import ARProcessor
from TOOLS_io import parse_par_file, parse_JSON, rsync, print_exe
from TOOLS_nullsar import fit_time_phase, fit_phase_offset, scale_freq_phase, plot_INIT, plot_OPT

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from injector import SCRIPT_inject_pulsars 
from injector.io_tools import merge_filterbanks, print_exe, FilterbankReader


class NullerProcess:
    def __init__(self, tag, processing_args, out_dir, work_dir, mode='INIT'):

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.processing_args_path = processing_args
        self.processing_args = parse_JSON(processing_args)['nullsar']

        self.tag = tag
        self.mode = mode
        self.processing_dir = f'{self.out_dir}/{self.tag}'

        self.ar_data = {}

    def injector_setup(self):
        self.check_SNR()

        self.get_data()
        self.extract_archive()
        self.create_injection_plan()

        self.gen_plots()

    def check_SNR(self):
        SNR_limit = self.processing_args.get('SNR_limit', 15)
        snr_path = f'{self.processing_dir}/02_INIT/INIT_fold_params.json'

        if self.mode == "INIT":
            self.SNR_record = {}
            for par_file  in list(self.processing_args['par_files']):
                psr_ID = Path(par_file).stem
                fits_path =  f'{self.processing_dir}/02_INIT/FOLDS/{psr_ID}_mode_INIT.fits'
                archive = ARProcessor(fits_path)
                SNR = archive.get_SNR()
                self.SNR_record[psr_ID] = {"SNR": SNR}

            with open(snr_path, 'w') as file:
                json.dump(self.SNR_record, file, indent=4)

        SNR_record = parse_JSON(snr_path)
        for par_file in list(self.processing_args['par_files']):
            psr_ID = Path(par_file).stem
            SNR = SNR_record[psr_ID]['SNR']
            if SNR < SNR_limit:
                self.processing_args['par_files'].remove(par_file)
                
        if len(self.processing_args['par_files']) == 0:
            print_exe('No pulsars to null.')
            sys.exit(0)

    def get_data(self):
        files_dir = f'{self.processing_dir}/01_FILES'

        if self.processing_args.get('filtool', False):
            data = glob.glob(f'{files_dir}/{self.tag}_FILTOOL*.fil')[0]
        else:
            with open(f'{files_dir}/files.txt') as f:
                data = [line.strip() for line in f if line.strip()]

        if type(data) == list:
            fb_names = [Path(fb).stem for fb in data]

            prefix = os.path.commonprefix(fb_names)
            self.new_fb_path = f"{self.work_dir}/{prefix}_MERGED.fil"

            for data_product in data:
                rsync(data_product, self.work_dir)
            
            new_data_paths = [f'{self.work_dir}/{Path(fb).name}' for fb in data]
            new_data_paths.sort()
            merge_filterbanks(new_data_paths, self.new_fb_path)

        else:
            prefix = Path(data).stem 
            self.new_fb_path = f"{self.work_dir}/{prefix}.fil"

            rsync(data, self.new_fb_path)

        ephem = self.processing_args['injection']['ephem']
        if ephem != 'builtin':
            rsync(ephem, self.work_dir)
            self.ephem = f'./{Path(ephem).name}'

    def extract_archive(self):
        par_files = self.processing_args['par_files']
        params_path = f'{self.processing_dir}/02_INIT/INIT_fold_params.json'
        
        for par_file in par_files:
            psr_ID = Path(par_file).stem
            if self.mode == 'INIT':
                fits_path = f'{self.processing_dir}/02_INIT/FOLDS/{psr_ID}_mode_INIT.fits'
                archive_INIT = ARProcessor(fits_path)
                out = self.parse_archive(psr_ID, archive_INIT, archive_INIT, params_path)

                save_path = f'{self.processing_dir}/02_INIT/MODELS/model_{psr_ID}.png'
                plot_INIT(save_path, archive_INIT, out)

            if self.mode == 'NULL':
                fits_path_INIT = f'{self.processing_dir}/02_INIT/FOLDS/{psr_ID}_mode_INIT.fits'
                fits_path_OPT = f'{self.processing_dir}/03_OPT/FOLDS/{psr_ID}_mode_OPTIMISE.fits'

                archive_INIT = ARProcessor(fits_path_INIT)
                archive_OPT = ARProcessor(fits_path_OPT)

                out = self.parse_archive(psr_ID, archive_INIT, archive_OPT, params_path)
                
                save_path = f'{self.processing_dir}/03_OPT/MODELS/OPT_{psr_ID}.png'
                plot_OPT(save_path, archive_INIT, archive_OPT, out)

        if self.mode == "INIT":
            with open(params_path, 'w') as file:
                for key, value in self.SNR_record.items():
                    if not self.ar_data.get(key, None):
                        self.ar_data[key] = value
                json.dump(self.ar_data, file, indent=4)

    def parse_archive(self, psr_ID, archive_INIT, archive_OPT, params_path):

        fb = FilterbankReader(self.new_fb_path, load_fb_stats=(128, 6))
        obs_len = fb.dt * fb.n_samples

        profile_path = f'{self.processing_dir}/02_INIT/MODELS/profile_{psr_ID}.npy'
        flux_time_path = f'{self.processing_dir}/02_INIT/MODELS/light_curve_{psr_ID}.npy'        
        if self.mode == 'INIT':
            SNR = archive_INIT.get_SNR()
            DM = archive_INIT.get_DM()
            freq_phase = archive_INIT.get_freq_phase()
            intensity_profile = archive_INIT.get_intensity_prof()
            
            time_phase = archive_INIT.get_time_phase()
            freq_deriv, phase_offset, phase_shift, SNR_fit, time, time_amp = fit_time_phase(time_phase, intensity_profile, obs_len)

            np.save(flux_time_path, np.stack([time, time_amp]))
            
            freq_phase_scaled = scale_freq_phase(freq_phase, intensity_profile)
            np.save(profile_path, freq_phase_scaled)

            self.ar_data[psr_ID] = {"SNR": SNR,  
                                    "DM": DM,
                                    "phase_offset": phase_offset-phase_shift, 
                                    "light_curve": flux_time_path,
                                    "profile": profile_path,
                                    "FX": freq_deriv}

            return freq_deriv, phase_offset, time, time_amp, obs_len

        elif self.mode == "NULL":
            init_ar_data =  parse_JSON(params_path)
            SNR = init_ar_data[psr_ID]['SNR']
            DM = init_ar_data[psr_ID]['DM']
            freq_deriv = init_ar_data[psr_ID]['FX']

            intensity_profile_INIT = archive_INIT.get_intensity_prof()
            intensity_profile_OPT = archive_OPT.get_intensity_prof()

            phase_offset, SNR_scale, fit_params = fit_phase_offset(intensity_profile_OPT, intensity_profile_INIT)

            self.ar_data[psr_ID] = {"SNR": SNR*SNR_scale,  
                                "DM": DM,
                                "phase_offset": phase_offset+init_ar_data[psr_ID]['phase_offset'], 
                                "light_curve": flux_time_path,
                                "profile": profile_path,
                                "FX": freq_deriv}

            return fit_params

        
            
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
                "phase_offset": float(self.ar_data[psr_ID]['phase_offset']),
                
                "P0_SNR": 1/float(params['F0']),
                "DM": self.ar_data[psr_ID]['DM'],
                "SNR": -self.ar_data[psr_ID]['SNR'],

                "profile": self.ar_data[psr_ID]['profile'],
                "light_curve": self.ar_data[psr_ID]['light_curve'],
                "polycos": par_file
            }

            for key, value in self.ar_data[psr_ID]['FX'].items():
                psr_dict[key] = value

            injection_plan['pulsars'].append(psr_dict)

        self.inject_file = f"{self.work_dir}/NULLSAR_inject_file_mode_{self.mode}.json"
        with open(self.inject_file, 'w') as file:
            json.dump(injection_plan, file, indent=4)

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
        report = f'{self.work_dir}/report_*.json'
        if not glob.glob(injected_fb):
            sys.exit(1)

        if self.mode == 'INIT':
            results_dir = f'{self.processing_dir}/02_INIT'
        elif self.mode == 'NULL':
            results_dir = f'{self.processing_dir}/03_OPT'


        rsync(injected_fb, results_dir)
        rsync(glob.glob(report)[0], results_dir)

        print_exe('Done.')


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Nullsar zapper',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='file tag')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--mode', metavar='int', type=str, required=False, default="init", help='folding mode (INIT, NULL)')
    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of cpus for injection')

    args = parser.parse_args()
    inj_exec = NullerProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.mode)
    inj_exec.injector_setup()
    inj_exec.run_injector(args.ncpus)
    inj_exec.transfer_products()