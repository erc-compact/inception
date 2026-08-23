import os
import glob
import argparse
import subprocess
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

        self.harmonic_log = {}

    def fold_setup(self):
        self.get_injection_report()
        self.transfer_data()
        self.create_zap_sting()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']
        self.cdm = self.injection_report['pulsars'][0]['cDM']

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]

        if self.processing_args['pulsarx_parfold_args'].get('transfer_TMP', True):
            inj_tools.rsync(data, self.work_dir)
            self.data = f'{self.work_dir}/{Path(data).name}'
        else:
            self.data = data

    def create_zap_sting(self):
        cmask = self.processing_args['pulsarx_parfold_args'].get('channel_mask', '')
        if cmask:
            cmask = cmask.strip()
            self.zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        else:
            self.zap_string = ''

    def set_blocksize(self, psr):
        default_plan = {'period_cutoff_s': [0.1, float('inf')], 'blocksize': [2, 10]}
        blocksize_plan = self.processing_args['pulsarx_parfold_args'].get('blocksize_plan', default_plan)

        period = psr['PX'][0]
        for period_max, block_size in zip(blocksize_plan['period_cutoff_s'], blocksize_plan['blocksize']):
            if period <= float(period_max):
                return block_size
            
    def set_nbins(self, psr):
        default_plan = {'period_cutoff_s': [0.1, float('inf')], 'nbins': [64, 128]}
        nbin_plan = self.processing_args['pulsarx_parfold_args'].get('nbin_plan', default_plan)

        period = psr['PX'][0]
        for period_max, nbins in zip(nbin_plan['period_cutoff_s'], nbin_plan['nbins']):
            if period <= float(period_max):
                return nbins
            
    def set_tsubints(self):
        fold_args = self.processing_args['pulsarx_parfold_args']
        obs_len = self.injection_report['injection_report']['obs_len']
        nsubint = fold_args.get('nsubint', 64)
        tsubint = fold_args['cmd'].get('tsubint', obs_len/64)
        if nsubint:
            tsubint = obs_len / nsubint
        return tsubint

    def get_parfile(self, psr):
        par_file =  f"{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr['ID']}.par"
        if glob.glob(par_file):
            return f'--parfile {par_file}'
        else:
            cand_file =  f"{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr['ID']}.candfile"
            return f'--candfile {cand_file}'

    def run_parfold(self, psr):
        fold_args = self.processing_args['pulsarx_parfold_args']

        psr_id = psr['ID']
        fold_file =  self.get_parfile(psr)
        block_size = self.set_blocksize(psr)
        nbins = self.set_nbins(psr)
        tsubint = self.set_tsubints()

        save_fits = '--saveimage' if fold_args.get('save_fits', False) else ''

        tmp_cwd = f'{self.work_dir}/process_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)

        template = f'{tmp_cwd}/TMP_template.template'
        Path(template).touch()

        cmd = f"{fold_args['mode']} -o {tmp_cwd}/{psr_id} --cdm {self.cdm} -f {self.data} --tsubint {tsubint} --output_width  --template {template} {fold_file} --blocksize {block_size} --nbin {nbins} {self.zap_string} {save_fits}"
    
        for flag in self.processing_args['pulsarx_parfold_args']['cmd_flags']:
            if flag in ['--output_width', '--saveimage']:
                pass
            cmd += f" {flag}"

        for key, value in self.processing_args['pulsarx_parfold_args']['cmd'].items():
            if key in ['nbin', 'blocksize', 'tsubint']:
                pass
            cmd += f" --{key} {value}"
        
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.png', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.ar', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.cands', self.work_dir)
        inj_tools.rsync(f'{tmp_cwd}/{psr_id}*.fits', self.work_dir)

        Path(template).unlink(missing_ok=True)

    def run_fold(self, ncpus):
        args = self.injection_report['pulsars']

        with Pool(ncpus) as p:
            p.map(self.run_parfold, args)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars'
        psr_ids = [arg['ID'] for arg in self.injection_report['pulsars']]

        if self.processing_args['pulsarx_parfold_args'].get('save_png', True):
            for pID in psr_ids:
                png = glob.glob(f'{self.work_dir}/{pID}*.png')
                if png:
                    os.rename(png[0], f"{Path(png[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.png")
            inj_tools.rsync(f'{self.work_dir}/*.png', results_dir)

        if self.processing_args['pulsarx_parfold_args'].get('save_ar', True):
            for pID in psr_ids:
                arc = glob.glob(f'{self.work_dir}/{pID}*.ar')
                if arc:
                    os.rename(arc[0], f"{Path(arc[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.ar")
            inj_tools.rsync(f'{self.work_dir}/*.ar', results_dir)

        if self.processing_args['pulsarx_parfold_args'].get('save_fits', False):
            for pID in psr_ids:
                fits = glob.glob(f'{self.work_dir}/{pID}*.fits')
                if fits:
                    os.rename(fits[0], f"{Path(fits[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.fits")
            inj_tools.rsync(f'{self.work_dir}/*.fits', results_dir)

        if self.processing_args['pulsarx_parfold_args'].get('save_cand', True):
            for pID in psr_ids:
                cands = glob.glob(f'{self.work_dir}/{pID}*.cands')
                if cands:
                    os.rename(cands[0], f"{Path(cands[0]).parent}/{pID}_{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}.cands")
            inj_tools.rsync(f'{self.work_dir}/*.cands', results_dir)

        if self.processing_args['pulsarx_parfold_args']['delete_inj_fb']:
            os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')

        if self.processing_args.get('pulsarx_candfold_args', {}).get('delete_inj_fb', False):
            check_cands = glob.glob(f'{self.out_dir}/inj_{self.injection_number:06}/inj_cands/*.png')
            if check_cands:
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