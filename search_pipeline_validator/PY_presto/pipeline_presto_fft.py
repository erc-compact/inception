import os
import glob
import argparse
import subprocess
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools


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

    def segment_root(self, dm):
        return f'{self.inj_id}_SEG_{self.seg_i}_{self.seg_n}_DS{self.downsample}_DM{dm:.2f}'

    def process_trial(self, dm):
        s_args = self.processing_args['presto_search_args']
        dat_dir = f'{self.results_dir}/processing/PRESTO/DAT'
        root = self.segment_root(dm)

        seg_dat = glob.glob(f'{dat_dir}/{root}.dat')
        seg_inf = glob.glob(f'{dat_dir}/{root}.inf')
        if not (seg_dat and seg_inf):
            inj_tools.print_exe(f'Missing dedispersed data for {root}, skipping.')
            return

        cwd = f'{self.work_dir}/DM{dm:.2f}'
        os.makedirs(cwd, exist_ok=True)
        inj_tools.rsync(seg_dat[0], cwd)
        inj_tools.rsync(seg_inf[0], cwd)

        cmd = f"realfft {cwd}/{root}.dat"
        cmd = inj_tools.add_cmd_args(cmd, s_args.get('realfft', {}))
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=cwd)

        if self.birdies:
            cmd = f"zapbirds -zap {self.birdies} {cwd}/{root}.fft"
            cmd = inj_tools.add_cmd_args(cmd, s_args.get('zapbirds', {}),
                                         skip_flags=['-zap'], skip_keys=['zapfile'])
            inj_tools.print_exe(cmd)
            subprocess.run(cmd, shell=True, cwd=cwd)

        if self.processing_args['presto_candfold_args'].get('fold_mode', 'filterbank') != 'dat':
            os.remove(f'{cwd}/{root}.dat')

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

        dat_dir = f'{self.results_dir}/processing/PRESTO/DAT'
        for dm in self.DM_list:
            for ext in ('.dat', '.inf'):
                f = f'{dat_dir}/{self.segment_root(dm)}{ext}'
                if os.path.exists(f):
                    os.remove(f)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto FFT (stage 2/3) for search pipeline validator',
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
