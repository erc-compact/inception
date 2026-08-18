import os
import sys
import glob
import argparse
import subprocess
from pathlib import Path

import search_pipeline_validator_TEST.pipeline_tools as inj_tools


class PulsarxFoldCandProcess:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

    def fold_setup(self):
        self.get_injection_report()
        self.transfer_data()
        self.create_zap_sting()
        self.get_candidates()

    def get_injection_report(self):
        report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")
        if data:
            inj_tools.rsync(data[0], self.work_dir)

            self.data = f'{self.work_dir}/{Path(data[0]).name}'
        else:
            sys.exit(0)

    def create_zap_sting(self):
        cmask = self.processing_args['pulsarx_candfold_args'].get('channel_mask', '')
        if cmask:
            cmask = cmask.strip()
            self.zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        else:
            self.zap_string = ''
    
    def get_candidates(self):
        processing_dir = f'{self.results_dir}/processing/PEASOUP'
        self.candfile = glob.glob(f"{processing_dir}/{self.process_tag}*.candfile")


    def set_tsubints(self):
        fold_args = self.processing_args['pulsarx_candfold_args']
        seg_div = self.process_tag.split('_')[-1]
        obs_len = self.injection_report['injection_report']['obs_len'] / seg_div

        nsubint = fold_args.get('nsubint', 64)
        tsubint = fold_args['cmd'].get('tsubint', obs_len/64)
        if nsubint:
            tsubint = obs_len / nsubint
        return tsubint


    def run_fold(self, ncpus):
        fold_args = self.processing_args['pulsarx_candfold_args']
        tsubint = self.set_tsubints()
        cmd = f"{fold_args['mode']} -t {ncpus} -o {self.work_dir}/ --tsubint {tsubint} -f {self.data} --template {fold_args['template']} --candfile {self.candfile} {self.zap_string}"
    
        for flag in self.processing_args['pulsarx_candfold_args']['cmd_flags']:
            cmd += f" {flag}"

        for key, value in self.processing_args['pulsarx_candfold_args']['cmd'].items():
            cmd += f" --{key} {value}"
        
        subprocess.run(cmd, shell=True)


    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_cands'
        if os.path.isfile(results_dir):
            os.remove(results_dir)
        os.makedirs(results_dir, exist_ok=True)

        for psr in self.injection_report['pulsars']:
            os.makedirs(f'{results_dir}/{psr}', exist_ok=True)

        if self.processing_args['pulsarx_candfold_args']['save_png']:
            pngs = glob.glob(f'{self.work_dir}/*.png')
            for png in pngs:
                psr_ID = png.split('_')[-1]
                os.rename(png, f"{Path(png).parent}/{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}_{psr_ID}")
            inj_tools.rsync(f'{self.work_dir}/*.png', results_dir)

        if self.processing_args['pulsarx_candfold_args']['save_ar']:
            arcs = glob.glob(f'{self.work_dir}/*.ar')
            for ar in arcs:
                psr_ID = ar.split('_')[-1]
                os.rename(ar, f"{Path(ar).parent}/{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}_{psr_ID}")
            inj_tools.rsync(f'{self.work_dir}/*.ar', results_dir)

        if self.processing_args['pulsarx_candfold_args']['save_cand']:
            inj_tools.rsync(f'{self.work_dir}/*.cands', results_dir)

        if self.processing_args['pulsarx_candfold_args']['delete_inj_fb']:
            check_par = glob.glob(f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/*.png')
            if check_par:
                os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Pulsarx candidate-folder for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='int', required=False, default=0, type=int, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = PulsarxFoldCandProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fold_exec.fold_setup()
    fold_exec.run_fold(args.ncpus)
    fold_exec.transfer_products()