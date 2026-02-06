import os
import csv
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
        df = pd.read_csv(candfile)
        all_cands = []
        for i, row in df.iterrows():
            all_cands.append([row['DM'], row['P(ms)'], row['candnum'], f"{presto_out_dir}/PRESTO_CANDS/{row['file']}.cand", 
                              row['SNR'], row['sigma'], row['numharm'], row['ipow'], row['cpow'], row['r'], row['z'], row['numhits']])

        self.cands = []
        for psr in self.injection_report['pulsars']:
            psr_DM = psr['DM']
            psr_p0 = psr['PX'][0]*1000
            for row in all_cands:
                if self.match(psr_DM, psr_p0, row):
                    self.cands.append((psr['ID'], *row))
                    break
                
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/inj_cands_PRESTO'
        os.makedirs(presto_out_dir, exist_ok=True)

        with open(f"{presto_out_dir}/candidates.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["ID", "DM", "P(ms)", "candnum", "file", 'SNR', 'sigma', 'numharm', 'ipow', 'cpow', 'r', 'z', 'numhits'])
            writer.writerows(self.cands)
        
    def match(self, psr_DM, psr_p0, row):
        DM, P0, *rest = row
        DM_tol = 0.2
        P0_ms_tol = 0.001
        if (abs(DM - psr_DM) <= DM_tol) and (abs(P0 - psr_p0) <= P0_ms_tol):
            return True
        else:
            return False

    def fold_candidate(self, candidates):
        ID, DM, P0, cand_i, candfile, *rest = candidates
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing'
        os.makedirs(presto_out_dir, exist_ok=True)

        f_args = self.processing_args['presto_search_args']

        if f_args.get('bary', False):
            bary = ""
        else:
            bary = "-topo"

        cwd = f"{self.work_dir}/{DM:.2f}_{cand_i}"
        out_file=f"{cwd}/{self.inj_id}_topo_DM{DM:.2f}"

        os.makedirs(cwd, exist_ok=True)
        cmd = f"prepfold {bary} -noxwin -n 64 -dm {DM:.2f} -o {out_file} -accelcand {cand_i} -accelfile {candfile} {self.mask} {self.data}"

        subprocess.run(cmd, shell=True, cwd=cwd)

        subprocess.run('ls', shell=True, cwd=cwd)

    def run_folds(self, ncpus):

        with Pool(ncpus) as p:
            p.map(self.fold_candidate, self.cands)


    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/inj_cands_PRESTO'
        os.makedirs(presto_out_dir, exist_ok=True)

        if self.processing_args['presto_fold_args']['save_png']:
            inj_tools.rsync(f'{self.work_dir}/*/*.png', presto_out_dir)

        if self.processing_args['presto_fold_args']['save_pfd']:
            inj_tools.rsync(f'{self.work_dir}/*/*.pfd', presto_out_dir)

        if self.processing_args['presto_fold_args']['delete_inj_fb']:
            check_par = glob.glob(f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/*.png')
            if check_par:
                os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')
            

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