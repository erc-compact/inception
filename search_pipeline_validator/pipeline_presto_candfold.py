import os
import glob
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool
import pandas as pd

import pipeline_tools as inj_tools


class PrestoFoldProcess:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def setup(self):
        self.get_injection_report()
        self.get_cands()
        self.transfer_data()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]
        inj_tools.rsync(data, self.work_dir)

        self.data = f'{self.work_dir}/{Path(data).name}'

    def get_cands(self):
        f_args = self.processing_args['presto_fold_args']
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing'

        if f_args['mask'] == 'rfifind':
            self.mask = f'-mask {presto_out_dir}/{self.inj_id}_rfifind.mask'
        elif f_args['mask']:
            self.mask = f"-mask {f_args['mask']}"
        else:
            self.mask = ''

        candfile = f"{presto_out_dir}/PRESTO_CANDS/PRESTO_candidates.txt"
        df = pd.readcsv(candfile)
        self.cands = []
        for row in df.iterrows():
            self.cands.append((row['DM'], row['candnum'], f"{presto_out_dir}/PRESTO_CANDS/{row['file']}.cand"))


    def fold_candidate(self, candidates):
        for (DM, cand_i, candfile) in candidates:
            results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
            presto_out_dir = f'{results_dir}/processing'
            os.makedirs(presto_out_dir, exist_ok=True)

            f_args = self.processing_args['presto_search_args']

            if f_args.get('bary', False):
                bary = ""
            else:
                bary = "-topo"

            out_file=f"{self.work_dir}/{self.inj_id}_topo_DM{DM:.2f}"
            cmd = f"prepfold {bary} -noxwin -dm {DM} -o {out_file} -accelcand {cand_i} -accelfile {candfile} {self.mask} {self.data}"

            subprocess.run(cmd, shell=True)

        subprocess.run("ls", shell=True)

    def run_folds(self, ncpus):

        with Pool(ncpus) as p:
            p.map(self.fold_candidate, self.cands)


    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing'
        out_file=f"{self.work_dir}/{self.inj_id}"

        if self.processing_args['presto_search_args']['save_png']:
            inj_tools.rsync(f'{out_file}*.png', presto_out_dir)

        if self.processing_args['presto_search_args']['save_pfd']:
            inj_tools.rsync(f'{out_file}*.pfd', presto_out_dir)
            

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='presto folder for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--threads', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = PrestoFoldProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fold_exec.setup()
    fold_exec.run_folds(args.threads)
    fold_exec.transfer_products()