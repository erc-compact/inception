import os
import glob
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

import pipeline_tools as inj_tools


class PrestoFoldParProcess:
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

        if self.processing_args['presto_parfold_args'].get('transfer_TMP', True):
            inj_tools.rsync(data, self.work_dir)
            self.data = f'{self.work_dir}/{Path(data).name}'
        else:
            self.data = data

    def set_nbins(self, psr):
        default_plan = {'period_cutoff_s': [0.1, float('inf')], 'nbins': [64, 128]}
        nbin_plan = self.processing_args['presto_parfold_args'].get('nbin_plan', default_plan)

        period = psr['PX'][0]
        for period_max, nbins in zip(nbin_plan['period_cutoff_s'], nbin_plan['nbins']):
            if period <= float(period_max):
                return nbins
            
    def get_mask(self):
        f_args = self.processing_args['presto_parfold_args']
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing/PRESTO'
        if f_args['mask'] == 'rfifind':
            return f'-mask {presto_out_dir}/{self.inj_id}_rfifind.mask'
        elif f_args['mask']:
            return f"-mask {f_args['mask']}"
        else:
            return ''

    def run_parfold(self, psr):
        fold_args = self.processing_args['presto_parfold_args']

        psr_id = psr['ID']
        par_file =  f"{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr['ID']}.par"
        if not glob.glob(par_file):
            return
        nbins = self.set_nbins(psr)
        mask = self.get_mask()
       
        tmp_cwd = f'{self.work_dir}/_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)

        cmd = f"prepfold -noxwin -o {tmp_cwd}/{psr_id} -filterbank {self.data} -par {par_file} -n {nbins} {mask}"
    
        for flag in self.processing_args['presto_parfold_args']['cmd_flags']:
            if flag in ['-noxwin']:
                pass
            cmd += f" {flag}"

        for key, value in self.processing_args['presto_parfold_args']['cmd'].items():
            if key in ['n']:
                pass
            cmd += f" -{key} {value}"
        
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*pfd*', self.work_dir)

    def run_fold(self, ncpus):
        args = self.injection_report['pulsars']

        with Pool(ncpus) as p:
            p.map(self.run_parfold, args)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars'
        psr_ids = [arg['ID'] for arg in self.injection_report['pulsars']]

        if self.processing_args['presto_parfold_args'].get('save_png', True):
            for pID in psr_ids:
                png = glob.glob(f'{self.work_dir}/{pID}*.png')
                if png:
                    os.rename(png[0], f"{Path(png[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}_pfd.png")
            inj_tools.rsync(f'{self.work_dir}/*.png', results_dir)

        if self.processing_args['presto_parfold_args'].get('save_pfd', True):
            for pID in psr_ids:
                pfd = glob.glob(f'{self.work_dir}/{pID}*.pfd')
                if pfd:
                    os.rename(pfd[0], f"{Path(pfd[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.pfd")
            inj_tools.rsync(f'{self.work_dir}/*.pfd', results_dir)

        if self.processing_args['presto_parfold_args'].get('save_bestprof', False):
            for pID in psr_ids:
                bestprof = glob.glob(f'{self.work_dir}/{pID}*.bestprof')
                if bestprof:
                    os.rename(bestprof[0], f"{Path(bestprof[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.bestprof")
            inj_tools.rsync(f'{self.work_dir}/*.bestprof', results_dir)

        if self.processing_args['presto_parfold_args'].get('save_polycos', True):
            for pID in psr_ids:
                polycos = glob.glob(f'{self.work_dir}/{pID}*.polycos')
                if polycos:
                    os.rename(polycos[0], f"{Path(polycos[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.polycos")
            inj_tools.rsync(f'{self.work_dir}/*.polycos', results_dir)

        if self.processing_args['presto_parfold_args'].get('save_ps', True):
            for pID in psr_ids:
                ps = glob.glob(f'{self.work_dir}/{pID}*.ps')
                if ps:
                    os.rename(ps[0], f"{Path(ps[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.ps")
            inj_tools.rsync(f'{self.work_dir}/*.ps', results_dir)

        if self.processing_args['presto_parfold_args']['delete_inj_fb']:
            os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')

        if self.processing_args.get('presto_candfold_args', {}).get('delete_inj_fb', False):
            check_cands = glob.glob(f'{self.out_dir}/inj_{self.injection_number:06}/inj_cands/PRESTO/*/*.png')
            if check_cands:
                os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto par-folder for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = PrestoFoldParProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fold_exec.fold_setup()
    fold_exec.run_fold(args.ncpus)
    fold_exec.transfer_products()