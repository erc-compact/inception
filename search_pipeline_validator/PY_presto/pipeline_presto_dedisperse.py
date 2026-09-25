import os
import glob
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from injector.io_tools import FilterbankReader


class PrestoDedisperseProcess:
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
        self.transfer_data()
        self.get_mask()
        self.get_segment()

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

    def transfer_data(self):
        data = glob.glob(f"{self.results_dir}/*_{self.inj_id}.fil")[0]
        if self.processing_args['presto_search_args'].get('transfer_TMP', True):
            inj_tools.rsync(data, self.work_dir)
            self.data = f'{self.work_dir}/{Path(data).name}'
        else:
            self.data = data

    def get_mask(self):
        s_args = self.processing_args['presto_search_args']
        presto_out_dir = f'{self.results_dir}/processing/PRESTO'
        if s_args['mask'] == 'rfifind':
            self.mask = f'-mask {presto_out_dir}/{self.inj_id}_rfifind.mask'
        elif s_args['mask']:
            self.mask = f"-mask {s_args['mask']}"
        else:
            self.mask = ''

    def get_segment(self):
        fb_reader = FilterbankReader(self.data, stats_samples=0)
        self.start_frac, self.numout = inj_tools.segment_samples(fb_reader.n_samples, self.seg_i,
                                                                 self.seg_n, self.downsample)

    def run_trial(self, dm):
        s_args = self.processing_args['presto_search_args']['prepdata']
        bary = '' if self.processing_args['presto_search_args'].get('bary', False) else '-nobary'

        cwd = f'{self.work_dir}/DM{dm:.2f}'
        os.makedirs(cwd, exist_ok=True)
        out_file = f'{cwd}/{self.inj_id}_SEG_{self.seg_i}_{self.seg_n}_DS{self.downsample}_DM{dm:.2f}'

        cmd = (f"prepdata {bary} -o {out_file} -dm {dm:.2f} -downsamp {self.downsample} "
               f"-start {self.start_frac} -numout {self.numout} {self.mask} {self.data}")

        cmd = inj_tools.add_cmd_args(cmd, s_args, skip_flags=['-nobary'],
                                     skip_keys=['o', 'dm', 'downsamp', 'start', 'numout', 'mask'])

        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=cwd)

    def run_dedisperse(self, ncpus):
        with Pool(ncpus) as p:
            p.map(self.run_trial, self.DM_list)

    def transfer_products(self):
        dat_dir = f'{self.results_dir}/processing/PRESTO/DAT'
        os.makedirs(dat_dir, exist_ok=True)

        inj_tools.rsync(f'{self.work_dir}/*/*.dat', dat_dir)
        inj_tools.rsync(f'{self.work_dir}/*/*.inf', dat_dir)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto dedisperser (stage 1/3) for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    dedisp_exec = PrestoDedisperseProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    dedisp_exec.setup()
    dedisp_exec.run_dedisperse(args.ncpus)
    dedisp_exec.transfer_products()
