import os
import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path
from multiprocessing import Pool

import pipeline_tools as inj_tools


class PulsarxFoldParProcess:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):
        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number
        self.rng = np.random.default_rng(self.injection_number)

    def fold_setup(self):
        self.get_injection_report()
        self.transfer_data()
        self.create_zap_sting()

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

    def create_zap_sting(self):
        cmask = self.processing_args['pulsarx_parfold_args'].get('channel_mask', '')
        if cmask:
            cmask = cmask.strip()
            self.zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        else:
            self.zap_string = ''

    def set_blocksize(self, psr):
        blocksize_plan = self.processing_args['pulsarx_parfold_args'].get('blocksize_plan', [2, 10, 0.1])
        input_blocksize = self.processing_args['pulsarx_parfold_args'].get('blocksize', False)
        if not input_blocksize:
            if psr['PX'][0] > blocksize_plan[2]:
                block_size = blocksize_plan[1]
            else:
                block_size = blocksize_plan[0]
        else:
            block_size = input_blocksize
        return block_size

    def get_parfile(self, psr):
        par_file =  f"{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr['ID']}.par"
        if glob.glob(par_file):
            harmonic_pars = self.processing_args['pulsarx_parfold_args'].get('harmonic_fold', None)
            if harmonic_pars:
                import pandas as pd
            
                max_duty_cycle = harmonic_pars.get('max_duty_cycle', 1)
                values = harmonic_pars.get('values', [1])
                weights =  harmonic_pars.get('weights', [1])

                if psr['duty_cycle'] <= max_duty_cycle:
                    
                    p = np.array(weights)/np.sum(weights)
                    harmonic = self.rng.choice(values, p=p)
                else:
                    harmonic = 1

                par_df = pd.read_csv(par_file, sep='\t', index_col=0, header=None)
                par_df.T.iloc[0]['F0'] = np.float64(par_df.T.iloc[0]['F0'])/harmonic

                new_par = f"{self.work_dir}/{psr['ID']}.par"
                par_df.to_csv(new_par, sep='\t', header=False)
                return new_par, '--parfile'
            else:
                return par_file, '--parfile' 
        else:
            cand_file =  f"{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr['ID']}.candfile"
            return cand_file, '--candfile'

    def run_parfold(self, psr):
        fold_args = self.processing_args['pulsarx_parfold_args']

        psr_id = psr['ID']
        par_file, mode_cmd =  self.get_parfile(psr)
        block_size = self.set_blocksize(psr)

        tmp_cwd = f'{self.work_dir}/process_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)
        cmd = f"{fold_args['mode']} -o {tmp_cwd}/{psr_id} -f {self.data} --template {fold_args['template']} {mode_cmd} {par_file} --blocksize {block_size} {self.zap_string}"
    
        for flag in self.processing_args['pulsarx_parfold_args']['cmd_flags']:
            cmd += f" {flag}"

        for key, value in self.processing_args['pulsarx_parfold_args']['cmd'].items():
            cmd += f" --{key} {value}"
        
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.png', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.ar', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.cands', self.work_dir)

    def run_fold(self, ncpus):
        args = self.injection_report['pulsars']

        with Pool(ncpus) as p:
            p.map(self.run_parfold, args)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars'

        if self.processing_args['pulsarx_parfold_args']['save_png']:
            inj_tools.rsync(f'{self.work_dir}/*.png', results_dir)
        if self.processing_args['pulsarx_parfold_args']['save_ar']:
            inj_tools.rsync(f'{self.work_dir}/*.ar', results_dir)
        if self.processing_args['pulsarx_parfold_args']['save_cand']:
            inj_tools.rsync(f'{self.work_dir}/*.cands', results_dir)
        if self.processing_args['pulsarx_parfold_args']['save_csv']:
            from candidate_tools import par_cand2csv

            output = f"{results_dir}/{self.processing_args['injection_args']['id']}_{self.inj_id}_parfold.csv"
            par_cand2csv(self.injection_report, self.work_dir, output)

        if self.processing_args['pulsarx_parfold_args']['delete_inj_fb']:
            os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Pulsarx par-folder for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = PulsarxFoldParProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fold_exec.fold_setup()
    fold_exec.run_fold(args.ncpus)
    fold_exec.transfer_products()