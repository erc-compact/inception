import os
import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path
from multiprocessing import Pool

import pipeline_tools as inj_tools


class PrestoSearchProcess:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def setup(self):
        self.get_injection_report()
        self.get_DDPlan()
        self.transfer_data()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def get_DDPlan(self):
        s_args = self.processing_args['presto_search_args']
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing'

        if s_args['mask'] == 'rfifind':
            self.mask = f'-mask {presto_out_dir}/{self.inj_id}_rfifind.mask'
        elif s_args['mask']:
            self.mask = f"-mask {s_args['mask']}"
        else:
            self.mask = ''

        if s_args['birdies'] == 'rfifind':
            self.birdies = f'-zapfile {presto_out_dir}/{self.inj_id}_birdies.txt'
        elif s_args['birdies']:
            self.birdies = f"-zapfile {s_args['birdies']}"
        else:
            self.birdies = ''

        ddplan_list = inj_tools.create_DDplan(s_args['ddplan'])
        DD_plan = []
        for i, d in enumerate(ddplan_list):
            n_trials = int(round((d.high_dm - d.low_dm) / d.dm_step))

            include_endpoint = (i == len(ddplan_list) - 1)
            if include_endpoint:
                n_trials += 1

            DM_values = np.linspace(d.low_dm, d.high_dm, n_trials, endpoint=include_endpoint)
            DD_plan.extend((dm, d.tscrunch) for dm in DM_values)

        self.DD_plan = DD_plan


    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]
        inj_tools.rsync(data, self.work_dir)

        self.data = f'{self.work_dir}/{Path(data).name}'

    def run_DM_trial(self, DM_trial):
        for (DM, down_sample) in DM_trial:
            results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
            presto_out_dir = f'{results_dir}/processing'
            os.makedirs(presto_out_dir, exist_ok=True)

            s_args = self.processing_args['presto_search_args']

            if s_args.get('bary', False):
                bary = ""
            else:
                bary = "-nobary"

            out_file=f"{self.work_dir}/{self.inj_id}_topo_DM{DM:.2f}"
            cmd=f"prepdata {bary} -o {out_file} -dm {DM} -downsamp {down_sample} {self.mask} {self.data}"
            subprocess.run(cmd, shell=True)

            cmd=f"realfft {out_file}.dat"
            subprocess.run(cmd, shell=True)

            cmd=f"zapbirds -zap {self.birdies} {out_file}.fft"
            subprocess.run(cmd, shell=True)
            
            if s_args.get('wmax', 0):
                wmax = f"-wmax {s_args['wmax']}" 
            else:
                wmax = ""

            cmd=f"accelsearch -numharm {s_args['numharm']} -zmax {s_args['zmax']} {wmax} {out_file}.fft"
            subprocess.run(cmd, shell=True)

            subprocess.run(f"rm -rf {out_file}.dat")
            subprocess.run(f"rm -rf {out_file}.fft")

        subprocess.run("ls", shell=True)

    def run_search(self, ncpus):

        with Pool(ncpus) as p:
            p.map(self.run_DM_trial, self.DD_plan)


    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing/PRESTO_CANDS'
        os.makedirs(presto_out_dir, exist_ok=True)

        inj_tools.rsync(f'{self.work_dir}/*.txt', presto_out_dir)
        inj_tools.rsync(f'{self.work_dir}/*.inf', presto_out_dir)
        inj_tools.rsync(f'{self.work_dir}/*.cand', presto_out_dir)
        inj_tools.rsync(f'{self.work_dir}/*ACCEL*', presto_out_dir)
            

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='presto search for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--threads', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    search_exec = PrestoSearchProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    search_exec.setup()
    search_exec.run_search(args.threads)
    search_exec.transfer_products()