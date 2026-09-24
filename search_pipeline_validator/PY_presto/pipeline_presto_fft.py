import os
import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from injector.io_tools import read_datfile, print_exe

from presto.infodata import infodata


class PrestoFFTProcess:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

    def setup(self):
        self.get_injection_report()
        self.parse_tag()
        self.get_trials()

    def get_injection_report(self):
        report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def parse_tag(self):
        self.batch_index = int(self.process_tag.split('_BATCH_')[-1])

    def get_trials(self):
        s_args = self.processing_args['presto_search_args']
        trials = inj_tools.build_dm_trials(s_args['ddplan'], s_args.get('inj_DM', True), self.injection_report)
        batches = inj_tools.batch_trials(trials, s_args.get('batch_size', 10))
        self.trials = batches[self.batch_index]

    def get_birdies(self):
        s_args = self.processing_args['presto_search_args']
        presto_out_dir = f'{self.results_dir}/processing/PRESTO'
        if s_args['birdies'] == 'rfifind':
            path = f'{presto_out_dir}/{self.inj_id}_birdies.txt'
            return f'-zapfile {path}' if os.path.exists(path) else ''
        elif s_args['birdies']:
            return f"-zapfile {s_args['birdies']}"
        else:
            return ''

    def process_trial(self, trial):
        dm, ds = trial
        s_args = self.processing_args['presto_search_args']
        dat_dir = f'{self.results_dir}/processing/PRESTO/DAT'

        full_root = f'{self.inj_id}_DS{ds}_DM{dm:.2f}'
        full_dat = glob.glob(f'{dat_dir}/{full_root}.dat')
        full_inf = glob.glob(f'{dat_dir}/{full_root}.inf')
        if not (full_dat and full_inf):
            print_exe(f'Missing dedispersed data for DM={dm:.2f} DS={ds}, skipping.')
            return

        cwd = f'{self.work_dir}/DS{ds}_DM{dm:.2f}'
        os.makedirs(cwd, exist_ok=True)
        inj_tools.rsync(full_dat[0], cwd)
        inj_tools.rsync(full_inf[0], cwd)

        base_info = infodata(f'{cwd}/{full_root}.inf')
        full_data = read_datfile(f'{cwd}/{full_root}.dat', nbits=32)
        n_total = len(full_data)

        fold_mode = self.processing_args['presto_candfold_args'].get('fold_mode', 'filterbank')
        birdies = self.get_birdies()

        for s_plan in s_args['segment_plan'].keys():
            n_seg = int(s_plan)
            for seg_i in range(n_seg):
                start = int(np.floor(seg_i * n_total / n_seg))
                end = int(np.floor((seg_i + 1) * n_total / n_seg))
                seg_data = full_data[start:end]

                seg_root = f'{self.inj_id}_SEG_{seg_i}_{n_seg}_DS{ds}_DM{dm:.2f}'
                seg_dat = f'{cwd}/{seg_root}.dat'
                seg_inf = f'{cwd}/{seg_root}.inf'

                seg_data.astype(np.float32).tofile(seg_dat)

                seg_info = infodata(f'{cwd}/{full_root}.inf')
                seg_info.N = len(seg_data)
                seg_info.epoch = base_info.epoch + (start * base_info.dt) / 86400.0
                seg_info.basenm = seg_root
                seg_info.to_file(seg_inf)

                cmd = f"realfft {seg_dat}"
                inj_tools.print_exe(cmd)
                subprocess.run(cmd, shell=True, cwd=cwd)

                if birdies:
                    cmd = f"zapbirds -zap {birdies} {cwd}/{seg_root}.fft"
                    inj_tools.print_exe(cmd)
                    subprocess.run(cmd, shell=True, cwd=cwd)

                if fold_mode != 'dat':
                    os.remove(seg_dat)

        os.remove(f'{cwd}/{full_root}.dat')

    def run_fft(self, ncpus):
        with Pool(ncpus) as p:
            p.map(self.process_trial, self.trials)

    def transfer_products(self):
        fft_dir = f'{self.results_dir}/processing/PRESTO/FFT'
        os.makedirs(fft_dir, exist_ok=True)
        inj_tools.rsync(f'{self.work_dir}/*/*.fft', fft_dir)
        inj_tools.rsync(f'{self.work_dir}/*/*_SEG_*.inf', fft_dir)

        if self.processing_args['presto_candfold_args'].get('fold_mode', 'filterbank') == 'dat':
            fold_dat_dir = f'{self.results_dir}/processing/PRESTO/FOLD_DAT'
            os.makedirs(fold_dat_dir, exist_ok=True)
            inj_tools.rsync(f'{self.work_dir}/*/*_SEG_*.dat', fold_dat_dir)

        # the transient, un-segmented per-trial .dat/.inf have now been consumed
        dat_dir = f'{self.results_dir}/processing/PRESTO/DAT'
        for dm, ds in self.trials:
            full_root = f'{self.inj_id}_DS{ds}_DM{dm:.2f}'
            for ext in ('.dat', '.inf'):
                f = f'{dat_dir}/{full_root}{ext}'
                if os.path.exists(f):
                    os.remove(f)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto segmenter/FFT (stage 2/3) for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='batch process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fft_exec = PrestoFFTProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fft_exec.setup()
    fft_exec.run_fft(args.ncpus)
    fft_exec.transfer_products()
