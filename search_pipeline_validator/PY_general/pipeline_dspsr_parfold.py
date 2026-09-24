import os
import glob
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

import pipeline_tools as inj_tools


class DSPSRFoldParProcess:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):
        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def fold_setup(self):
        self.get_injection_report()
        self.transfer_data()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']
        self.cdm = self.injection_report['pulsars'][0]['cDM']

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]

        if self.processing_args['dspsr_parfold_args'].get('transfer_TMP', True):
            inj_tools.rsync(data, self.work_dir)
            self.data = f'{self.work_dir}/{Path(data).name}'
        else:
            self.data = data

    def set_nbins(self, psr):
        default_plan = {'period_cutoff_s': [0.1, float('inf')], 'nbins': [64, 128]}
        nbin_plan = self.processing_args['dspsr_parfold_args'].get('nbin_plan', default_plan)

        period = psr['PX'][0]
        for period_max, nbins in zip(nbin_plan['period_cutoff_s'], nbin_plan['nbins']):
            if period <= float(period_max):
                return nbins
    
    def set_tsubints(self):
        fold_args = self.processing_args['dspsr_parfold_args']
        obs_len = self.injection_report['injection_report']['obs_len']
        nsubint = fold_args.get('nsubint', 64)
        tsubint = obs_len / nsubint
        return tsubint
    
    def run_DSPSR(self, psr, par_file, tmp_cwd):
        dspsr_args = self.processing_args['dspsr_parfold_args']['dspsr']

        nbins = self.set_nbins(psr)
        tsubint = self.set_tsubints()

        psr_id = psr['ID']
        cmd = f"dspsr -A -O {tmp_cwd}/{psr_id} -E {par_file} -b {nbins} -L {tsubint}"
    
        for flag in dspsr_args['cmd_flags']:
            if flag in ['-A']:
                pass
            cmd += f" {flag}"

        for key, value in dspsr_args['cmd'].items():
            if key in ['b', 'L']:
                pass
            cmd += f" -{key} {value}"
        
        cmd += f' {self.data}'
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

    def run_PDMP(self, psr, tmp_cwd):
        pdmp_args = self.processing_args['dspsr_parfold_args']['pdmp']

        psr_id = psr['ID']
        cmd = f"pdmp -g {tmp_cwd}/{psr_id}.png/PNG"
    
        for flag in pdmp_args['cmd_flags']:
            if flag in ["-v"]:
                pass
            cmd += f" {flag}"

        for key, value in pdmp_args['cmd'].items():
            if key in []:
                pass
            cmd += f" -{key} {value}"
        
        cmd += f' {tmp_cwd}/{psr_id}.ar'
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

        inj_tools.rsync(f'{tmp_cwd}/{psr_id}.ar', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/{psr_id}.png', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/pdmp.best', f'{self.work_dir}/{psr_id}.best')
            
    def run_parfold(self, psr):
        psr_id = psr['ID']
        par_file =  f"{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr['ID']}.par"
        if not glob.glob(par_file):
            return
        
        tmp_cwd = f'{self.work_dir}/_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)

        self.run_DSPSR(psr, par_file, tmp_cwd)

        self.run_PDMP(psr, tmp_cwd)

    def run_fold(self, ncpus):
        args = self.injection_report['pulsars']

        with Pool(ncpus) as p:
            p.map(self.run_parfold, args)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars'
        psr_ids = [arg['ID'] for arg in self.injection_report['pulsars']]

        if self.processing_args['dspsr_parfold_args'].get('save_png', True):
            for pID in psr_ids:
                png = glob.glob(f'{self.work_dir}/{pID}*.png')
                if png:
                    os.rename(png[0], f"{Path(png[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}_dspsr.png")
            inj_tools.rsync(f'{self.work_dir}/*.png', results_dir)

        if self.processing_args['dspsr_parfold_args'].get('save_ar', True):
            for pID in psr_ids:
                ar = glob.glob(f'{self.work_dir}/{pID}*.ar')
                if ar:
                    os.rename(ar[0], f"{Path(ar[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}_dspsr.ar")
            inj_tools.rsync(f'{self.work_dir}/*.ar', results_dir)

        if self.processing_args['dspsr_parfold_args'].get('save_best', True):
            for pID in psr_ids:
                file = glob.glob(f'{self.work_dir}/{pID}.best')
                if file:
                    os.rename(file[0], f"{Path(file[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}_dspsr.best")
            inj_tools.rsync(f'{self.work_dir}/*.best', results_dir)

        if self.processing_args['dspsr_parfold_args']['delete_inj_fb']:
            os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='DSPSR par-folder for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = DSPSRFoldParProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fold_exec.fold_setup()
    fold_exec.run_fold(args.ncpus)
    fold_exec.transfer_products()