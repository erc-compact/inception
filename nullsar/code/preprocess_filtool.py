import os
import glob
import argparse
import subprocess
from pathlib import Path

from nullsar_tools import parse_JSON, rsync, print_exe


class FiltoolProcess:
    def __init__(self, tag, processing_args, out_dir, work_dir):

        self.processing_args_path = processing_args
        self.processing_args = parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.tag = tag

    def transfer_data(self):
        files_dir = f'{self.out_dir}/PROCESSING/{self.tag}/01_FILES'

        with open(f'{files_dir}/files.txt') as f:
            self.fb = [line.strip() for line in f if line.strip()]
        
        for fb_file in self.fb:
            rsync(fb_file, self.work_dir)

        self.data = " ".join([f'{self.work_dir}/{Path(data).name}' for data in self.fb])

    def run_filtool(self, threads):
        rootname = f"{self.tag}_FILTOOL" 

        cmd = f"filtool -t {threads} -o {self.work_dir}/{rootname} -f {self.data}"

        for flag in self.processing_args['filtool_args']['cmd_flags']:
            cmd += f" {flag}"

        for key, value in self.processing_args['filtool_args']['cmd'].items():
            cmd += f" --{key} {value}"
        
        subprocess.run(cmd, shell=True)

    def transfer_products(self):
        results_dir = f'{self.out_dir}/PROCESSING/{self.tag}/01_FILES'

        data_product = glob.glob(f'{self.work_dir}/*_FILTOOL*.fil')[0]
        rsync(data_product, results_dir)

    def process(self, threads):

        check_fb = glob.glob(f'{self.out_dir}/PROCESSING/{self.tag}/01_FILES/*_FILTOOL*.fil')
        if check_fb:
            print_exe('Filterbank already FILTOOLED.')
        else:
            self.transfer_data()
            self.run_filtool(threads)
            self.transfer_products()
            

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='NULLSAR - preprocessing filtool',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='file tag')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--threads', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    fil_exec = FiltoolProcess(args.tag, args.processing_args, args.out_dir, args.work_dir)

    fil_exec.process(args.threads)