import re
import sys
import json
import argparse
import subprocess
import pandas as pd
from pathlib import Path


class CandExec:
    def __init__(self, candidate_files, search_args, injection_report, output_dir):
        self.out = output_dir
        self.candidate_files = candidate_files
        args = self.parse_JSON(search_args)
        self.inj_report = self.parse_JSON(injection_report)
        self.multi_beam, self.ID = args['multi_beam_args'], args['processing_id']

    def parse_JSON(self, json_file):
        try:
            with open(json_file, 'r') as file:
                pars = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {json_file}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {json_file} using JSON.')
        else:
            return pars
        
    def add_beam_info(self):
        def get_beam_name(file_path):
            path_obj = Path(file_path)
            if path_obj.stem.split('_')[0] == 'data':
                return re.search(r'[ci]fbf\d{5}', self.inj_report['injection']['fb']).group()
            else:
                return path_obj.parent.stem 

        fold_candidates = pd.read_csv(f'{self.out}/_good_cands_to_fold.csv')
        fold_candidates['beam_id'] = [get_beam_name(file_path) for file_path in fold_candidates['file']]
        
        fold_candidates.to_csv(f'{self.out}/_good_cands_to_fold_with_beam.csv')

    def run_cmd(self):
        cmd = f"candidate_filter.py -i {self.candidate_files} -o {self.out}/ --threshold {self.multi_beam['snr_cutoff']} " \
              f"--p_tol {self.multi_beam['p_tol']} --dm_tol {self.multi_beam['dm_tol']} " \
              "-c /home/psr/software/candidate_filter/candidate_filter/default_config.json " \
              "--rfi /home/psr/software/candidate_filter/candidate_filter/known_rfi.txt"

        subprocess.run(cmd, shell=True)

        self.add_beam_info()

    

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate filter for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--candidates_file', metavar='file', required=True, help='file containing peasoup xml names')
    parser.add_argument('--injection_report', metavar='file', required=True, help='JSON file with inject pulsar records')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    cand_exec = CandExec(args.candidates_file, args.search_args, args.injection_report, args.output)
    cand_exec.run_cmd()

