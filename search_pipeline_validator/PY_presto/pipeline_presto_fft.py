import os
import glob
import argparse
import subprocess
import numpy as np
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools

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
        self.dat_dir = f'{self.results_dir}/processing/PRESTO/DAT'

    def setup(self):
        self.get_injection_report()
        self.parse_tag()
        self.get_DM_list()
        self.get_birdies()

    def get_injection_report(self):
        report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def parse_tag(self):
        self.downsample, self.seg_i, self.seg_n = inj_tools.parse_process_tag(self.process_tag)

    def get_DM_list(self):
        s_args = self.processing_args['presto_search_args']
        self.DM_list = inj_tools.build_dm_list(s_args['ddplan'], self.downsample,
                                               s_args.get('inj_DM', True), self.injection_report)

    def get_birdies(self):
        s_args = self.processing_args['presto_search_args']
        presto_out_dir = f'{self.results_dir}/processing/PRESTO'
        if s_args['birdies'] == 'rfifind':
            path = f'{presto_out_dir}/{self.inj_id}_birdies.txt'
            self.birdies = f'-zapfile {path}' if os.path.exists(path) else ''
        elif s_args['birdies']:
            self.birdies = f"-zapfile {s_args['birdies']}"
        else:
            self.birdies = ''

    def full_root(self, dm):
        return f'{self.inj_id}_DS{self.downsample}_DM{dm:.2f}'

    def segment_root(self, dm):
        return f'{self.inj_id}_SEG_{self.seg_i}_{self.seg_n}_DS{self.downsample}_DM{dm:.2f}'

    def cut_segment(self, dm, cwd):
        """Read only this segment's slice out of the full .dat (the file is shared
        with every other segment of this downsample, so it is never copied)."""
        full_root = self.full_root(dm)
        full_dat = f'{self.dat_dir}/{full_root}.dat'
        full_inf = f'{self.dat_dir}/{full_root}.inf'

        if not (os.path.exists(full_dat) and os.path.exists(full_inf)):
            inj_tools.print_exe(f'Missing dedispersed data for {full_root}, skipping.')
            return None

        n_total = os.path.getsize(full_dat) // 4
        if n_total == 0:
            inj_tools.print_exe(f'Dedispersed data for {full_root} is empty, skipping.')
            return None

        start, count = inj_tools.segment_samples(n_total, self.seg_i, self.seg_n)

        seg_root = self.segment_root(dm)
        with open(full_dat, 'rb') as f:
            f.seek(start * 4)
            data = np.fromfile(f, dtype=np.float32, count=count)

        if data.size == 0:
            inj_tools.print_exe(f'Empty segment for {seg_root}, skipping.')
            return None

        data.tofile(f'{cwd}/{seg_root}.dat')

        seg_info = infodata(full_inf)
        seg_info.N = data.size
        seg_info.epoch = seg_info.epoch + (start * seg_info.dt) / 86400.0
        seg_info.basenm = seg_root
        seg_info.to_file(f'{cwd}/{seg_root}.inf')

        return seg_root

    def process_trial(self, dm):
        s_args = self.processing_args['presto_search_args']

        cwd = f'{self.work_dir}/DM{dm:.2f}'
        os.makedirs(cwd, exist_ok=True)

        seg_root = self.cut_segment(dm, cwd)
        if seg_root is None:
            return

        cmd = f"realfft {cwd}/{seg_root}.dat"
        cmd = inj_tools.add_cmd_args(cmd, s_args.get('realfft', {}), skip_flags=['-delete'])
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=cwd)

        fft_file = f'{cwd}/{seg_root}.fft'
        if (not os.path.exists(fft_file)) or os.path.getsize(fft_file) == 0:
            raise RuntimeError(f'realfft produced no output for {seg_root}.')

        if self.birdies:
            cmd = f"zapbirds -zap {self.birdies} {fft_file}"
            cmd = inj_tools.add_cmd_args(cmd, s_args.get('zapbirds', {}),
                                         skip_flags=['-zap'], skip_keys=['zapfile'])
            inj_tools.print_exe(cmd)
            subprocess.run(cmd, shell=True, cwd=cwd)

        if self.processing_args['presto_candfold_args'].get('fold_mode', 'filterbank') != 'dat':
            os.remove(f'{cwd}/{seg_root}.dat')

    def run_fft(self, ncpus):
        with Pool(ncpus) as p:
            p.map(self.process_trial, self.DM_list)

    def transfer_products(self):
        fft_dir = f'{self.results_dir}/processing/PRESTO/FFT'
        os.makedirs(fft_dir, exist_ok=True)
        inj_tools.rsync(f'{self.work_dir}/*/*.fft', fft_dir)
        inj_tools.rsync(f'{self.work_dir}/*/*.inf', fft_dir)

        if self.processing_args['presto_candfold_args'].get('fold_mode', 'filterbank') == 'dat':
            fold_dat_dir = f'{self.results_dir}/processing/PRESTO/FOLD_DAT'
            os.makedirs(fold_dat_dir, exist_ok=True)
            inj_tools.rsync(f'{self.work_dir}/*/*.dat', fold_dat_dir)
            inj_tools.rsync(f'{self.work_dir}/*/*.inf', fold_dat_dir)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto segmenter/FFT (stage 2/3) for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fft_exec = PrestoFFTProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fft_exec.setup()
    fft_exec.run_fft(args.ncpus)
    fft_exec.transfer_products()
