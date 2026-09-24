import os
import glob
import argparse
import subprocess

import pipeline_tools as inj_tools


class Classifier:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        c_args = self.processing_args['classifier_args']
        self.tag = c_args['tag']
        self.python = c_args['python']
        self.script = c_args['script']
        self.model_dir = c_args['model_dir']
        self.file_type = c_args['file_type']

    def classify(self):
        attempts = 0
        
        if self.file_type == 'pulsarx':
            cand_dir = f'{self.results_dir}/inj_cands/PEASOUP/{self.process_tag}' 
            extension = 'ar'
        elif self.file_type == 'presto':
            cand_dir = f'{self.results_dir}/inj_cands/PRESTO/{self.process_tag}' 
            extension = 'pfd'

        ar_files = glob.glob(f'{cand_dir}/*.{extension}')
        results_file = f'{cand_dir}/CLASSIFIED_{self.tag}.csv'

        if len(ar_files) != 0:
            while (not os.path.exists(results_file)) and (attempts != 3):
                attempts += 1
                cmd = f"{self.python} {self.script} --cand_dir={cand_dir} --model_dir={self.model_dir} --extension={extension} --output_file={results_file}"
                subprocess.run(cmd, shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='classifier for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    cl_exec = Classifier(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    cl_exec.classify()