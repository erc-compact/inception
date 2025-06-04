import os
import glob
import json
import argparse
import subprocess
from pathlib import Path

import pipeline_tools as inj_tools


class FiltoolProcess:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def filtool_setup(self):
        self.get_injection_report()
        self.transfer_data()
        self.create_filplan()

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

    def create_filplan(self):
        filtool_args = self.processing_args['filtool_args']
        filplan_args = filtool_args['filplan']
        filplan_file = f'{self.work_dir}/filplan.json'

        with open(filplan_file, 'w') as filplan:
            plans = []
            for tscrunch in filtool_args['tscrunch']:
                plans.append({'time_downsample': int(tscrunch),
                              'frequency_downsample': filplan_args['frequency_downsample'],
                              'baseline_width': filplan_args['baseline_width'],
                              'dataout_mean': filplan_args['outmean'],
                              'dataout_std': filplan_args['outstd'],
                              'dataout_nbits': filplan_args['nbits'],
                              'rfi_flags': filplan_args.get('rfi_flags', '')})
            json.dump(plans, filplan)

        self.filplan_file = filplan_file

    def run_filtool(self, threads):
        rootname = f"{Path(self.data).stem}_FILTOOL" 

        cmd = f"filtool -t {threads} --filplan {self.filplan_file}  -o {self.work_dir}/{rootname} -f {self.data}"

        for flag in self.processing_args['filtool_args']['cmd_flags']:
            cmd += f" {flag}"

        for key, value in self.processing_args['filtool_args']['cmd'].items():
            cmd += f" --{key} {value}"
        
        subprocess.run(cmd, shell=True)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        filtool_out_dir = f'{results_dir}/processing'
        os.makedirs(filtool_out_dir, exist_ok=True)

        if self.processing_args['filtool_args']['save_filtool_fb']:
            data_products = glob.glob(f'{self.work_dir}/*_FILTOOL*.fil')
            for filename in data_products:
                inj_tools.rsync(filename, filtool_out_dir)

        if self.processing_args['filtool_args']['save_filplan']:
            inj_tools.rsync(self.filplan_file, filtool_out_dir)

        if self.processing_args['filtool_args']['delete_inj_fb']:
            os.remove(f'{self.out_dir}/inj_{self.injection_number:06}/{Path(self.data).name}')
            

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='filtool for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--threads', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    fil_exec = FiltoolProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fil_exec.filtool_setup()
    fil_exec.run_filtool(args.threads)
    fil_exec.transfer_products()