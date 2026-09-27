import os
import glob
import argparse
from pathlib import Path

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools


class PrestoCleanup:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):
        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

    def clean(self):
        s_args = self.processing_args['presto_search_args']
        presto_dir = f'{self.results_dir}/processing/PRESTO'

        files_dir = f'{presto_dir}/FILES'

        if not s_args.get('save_dat', False):
            full_obs = [f for f in glob.glob(f'{files_dir}/*') if '_SEG_' not in Path(f).name]
            freed = sum(os.path.getsize(f) for f in full_obs)
            for f in full_obs:
                os.remove(f)
            inj_tools.print_exe(f'Removed dedispersed time series: {freed/1e9:.2f} GB.')

        if os.path.isdir(files_dir) and not os.listdir(files_dir):
            os.rmdir(files_dir)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto intermediate-product cleanup for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    clean_exec = PrestoCleanup(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    clean_exec.clean()
