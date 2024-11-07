import os
import sys
import json
import argparse
import subprocess

class FoldExec:
    def __init__(self, fb, search_args, output):
        self.out = output
        self.fb = fb
        args = self.parse_JSON(search_args)
        self.fold_args = args['fold_args']

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
        TEMPLATE = "/home/psr/software/PulsarX/include/template/meerkat_fold.template"

        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v -t 12 --candfile {} -n {} {} {} --template {TEMPLATE}" \
                f"--clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand -f {self.fb} --rfi zdot {}" \
                f"--fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']}"

        subprocess.run(cmd, shell=True)
        """
        pred_file, nsubband, nbins_string, beam_tag, TEMPLATE, subint_length, input_filenames, zap_string, fscrunch, tscrunch)
        """



if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate filter for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--candidates_file', metavar='file', required=True, help='file containing peasoup candidates ')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    fold_exec = FoldExec(args.output)
    fold_exec.run_cmd()