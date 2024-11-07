import os
import sys
import json
import argparse
import subprocess

class FoldExec:
    def __init__(self, output):
        self.out = output



    def run_cmd(self):

        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v -t 12 --candfile {} -n {} {} {} --template {} --clfd 8 -L {} --fillPatch rand -f {} --rfi zdot {} --fd {} --td {}"

        subprocess.run(cmd)
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