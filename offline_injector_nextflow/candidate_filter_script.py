import os
import sys
import json
import argparse
import subprocess


class CandExec:
    def __init__(self, candidate_files, search_args, output_dir):
        self.out = output_dir
        self.candidate_files = candidate_files
        args = self.parse_JSON(search_args)
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
        
    def run_cmd(self):

        cmd = f"candidate_filter.py -i {self.candidate_files} -o {self.out}/{self.ID} --threshold {self.multi_beam['snr_cutoff']}" \
              f"--p_tol {self.multi_beam['p_tol']} --dm_tol {self.multi_beam['dm_tol']} " \
              "-c /home/psr/software/candidate_filter/candidate_filter/default_config.json" \
              "--rfi /home/psr/software/candidate_filter/candidate_filter/known_rfi.txt"

        subprocess.run(cmd, shell=True)
    

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate filter for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--fb', metavar='file', required=True, help='Injected filterbank file')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    cand_exec = CandExec(args.fb, args.search_args, args.output)
    cand_exec.run_cmd()
