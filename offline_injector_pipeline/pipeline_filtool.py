import os
import json
import argparse
import subprocess
from pathlib import Path

from pipeline_tools import PipelineTools


class FiltoolExec(PipelineTools):
    def __init__(self, injection_number, search_args, out_dir, ncpus):
        super().__init__(search_args)
        self.work_dir = os.getcwd()
        self.out_dir = out_dir
        self.injection_number = injection_number
        self.fb = self.get_fb()
        self.ncpus = ncpus

        self.define_defaults()

    def get_fb(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        for filename in os.listdir(results_dir):
            if f'_{self.data_ID}_' in filename:
                subprocess.run(f"rsync -Pav {results_dir}/{filename} {self.work_dir}", shell=True)
                return f'{self.work_dir}/{Path(filename).name}'

    def define_defaults(self):
        self.pars = {
            'zapping_threshold': 4.0, 
            'segment_length': 2.0, 
            'baseline': '0 0', 
            'fillPatch': 'rand',
            'nbits': 8,
            'outmean': 128.0,
            'outstd': 6.0,
            'frequency_downsample': 1,
            'baseline_width': 0.0
        }
    
    def create_filplan(self, ddplan):
        filplan_file = os.path.join(self.work_dir, "filplan.json")

        with open(filplan_file, 'w') as filplan:
            plans = []
            for dm_range in ddplan:
                plans.append({"time_downsample": int(dm_range.tscrunch),
                            "frequency_downsample": self.pars['frequency_downsample'],
                            "baseline_width": self.pars['baseline_width'],
                            "dataout_mean": self.pars['outmean'],
                            "dataout_std": self.pars['outstd'],
                            "dataout_nbits": self.pars['nbits'],
                            "rfi_flags": ""})
            json.dump(plans, filplan)
        return filplan_file

    def run_cmd(self):
        ddplan = self.create_DDplan()
        filplan_file = self.create_filplan(ddplan)

        fscrunch = self.processing_args.get('fscrunch', 1) 
        rfi_flags = self.processing_args.get('rfi_flags', 'zdot')
        
        fb_root = Path(self.fb).stem.split('_')
        rootname = f"downsampled_{fb_root[-2]}_{fb_root[-1]}" 

        cmd = f"filtool -v -t {self.ncpus} --zapthre {self.pars['zapping_threshold']} --baseline {self.pars['baseline']} -l {self.pars['segment_length']} \
                 --filplan {filplan_file} --fd {fscrunch} --fillPatch {self.pars['fillPatch']} -z {rfi_flags} -o {self.work_dir}/{rootname} -f {self.fb}"
        
        subprocess.run(cmd, shell=True)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        filtool_out_dir = f'{results_dir}/processing'
        os.mkdir(filtool_out_dir)

        for filename in os.listdir(self.work_dir):
            if 'downsampled' in filename:
                subprocess.run(f"rsync -Pav {filename} {filtool_out_dir}", shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='filtool for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    parser.add_argument('--ncpus', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    fil_exec = FiltoolExec(args.injection_number, args.search_args, args.out_dir, args.ncpus)
    fil_exec.run_cmd()
    fil_exec.transfer_products()