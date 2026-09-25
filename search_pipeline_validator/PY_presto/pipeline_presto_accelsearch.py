import os
import glob
import shutil
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools


class PrestoAccelsearchProcess:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        self.fft_dir = f'{self.results_dir}/processing/PRESTO/FFT'

    def setup(self):
        self.get_injection_report()
        self.parse_tag()
        self.get_DM_list()
        self.get_jobs()

    def get_injection_report(self):
        report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def parse_tag(self):
        self.downsample, self.seg_i, self.seg_n = inj_tools.parse_process_tag(self.process_tag)
        self.seg_args = self.processing_args['presto_search_args']['segment_plan'][f'{self.seg_n}']

    def get_DM_list(self):
        s_args = self.processing_args['presto_search_args']
        self.DM_list = inj_tools.build_dm_list(s_args['ddplan'], self.downsample,
                                               s_args.get('inj_DM', True), self.injection_report)

    def segment_root(self, dm):
        return f'{self.inj_id}_SEG_{self.seg_i}_{self.seg_n}_DS{self.downsample}_DM{dm:.2f}'

    def get_jobs(self):
        # accelsearch reads the .fft and writes its ACCEL products alongside it,
        # so it runs against FFT/ directly rather than staging a copy of every
        # .fft (which is the same size as the time series it came from)
        self.jobs = []
        for dm in self.DM_list:
            root = self.segment_root(dm)
            if os.path.exists(f'{self.fft_dir}/{root}.fft') and os.path.exists(f'{self.fft_dir}/{root}.inf'):
                self.jobs.append(root)
            else:
                inj_tools.print_exe(f'No FFT found for {root}, skipping.')

    def run_accelsearch(self, root):
        s_args = self.processing_args['presto_search_args']

        wmax = f"-wmax {self.seg_args['wmax']}" if self.seg_args.get('wmax', 0) else ''
        sigma = f"-sigma {self.seg_args['sigma']}" if self.seg_args.get('sigma', None) else ''

        cmd = (f"accelsearch -numharm {self.seg_args['numharm']} -zmax {self.seg_args['zmax']} "
               f"{wmax} {sigma} {self.fft_dir}/{root}.fft")

        cmd = inj_tools.add_cmd_args(cmd, s_args.get('accelsearch', {}),
                                     skip_keys=['numharm', 'zmax', 'wmax', 'sigma'])

        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=self.work_dir)

    def run_search(self, ncpus):
        with Pool(ncpus) as p:
            p.map(self.run_accelsearch, self.jobs)

    def transfer_products(self):
        accel_dir = f'{self.results_dir}/processing/PRESTO/ACCEL'
        os.makedirs(accel_dir, exist_ok=True)

        save_fft = self.processing_args['presto_search_args'].get('save_fft', False)

        for root in self.jobs:
            # ACCEL products land next to the .fft; move only this job's files so
            # concurrent segments never touch each other's outputs
            for product in glob.glob(f'{self.fft_dir}/{root}_ACCEL_*'):
                shutil.move(product, f'{accel_dir}/{Path(product).name}')

            inf = f'{self.fft_dir}/{root}.inf'
            if os.path.exists(inf):
                shutil.copy(inf, f'{accel_dir}/{root}.inf')

            if not save_fft:
                for ext in ('.fft', '.inf'):
                    f = f'{self.fft_dir}/{root}{ext}'
                    if os.path.exists(f):
                        os.remove(f)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto accelsearch (stage 3/3) for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    search_exec = PrestoAccelsearchProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    search_exec.setup()
    search_exec.run_search(args.ncpus)
    search_exec.transfer_products()
