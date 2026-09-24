import os
import glob
import argparse
import subprocess
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools


class PrestoAccelsearchProcess:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

    def setup(self):
        self.get_injection_report()
        self.parse_tag()
        self.get_trials()
        self.transfer_data()

    def get_injection_report(self):
        report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def parse_tag(self):
        self.batch_index = int(self.process_tag.split('_BATCH_')[-1])

    def get_trials(self):
        s_args = self.processing_args['presto_search_args']
        trials = inj_tools.build_dm_trials(s_args['ddplan'], s_args.get('inj_DM', True), self.injection_report)
        batches = inj_tools.batch_trials(trials, s_args.get('batch_size', 10))
        self.trials = batches[self.batch_index]

    def transfer_data(self):
        fft_dir = f'{self.results_dir}/processing/PRESTO/FFT'
        segment_plan = self.processing_args['presto_search_args']['segment_plan']

        self.jobs = []
        for dm, ds in self.trials:
            for s_plan, seg_args in segment_plan.items():
                n_seg = int(s_plan)
                for seg_i in range(n_seg):
                    root = f'{self.inj_id}_SEG_{seg_i}_{n_seg}_DS{ds}_DM{dm:.2f}'
                    fft_file = glob.glob(f'{fft_dir}/{root}.fft')
                    inf_file = glob.glob(f'{fft_dir}/{root}.inf')
                    if fft_file and inf_file:
                        cwd = f'{self.work_dir}/{root}'
                        os.makedirs(cwd, exist_ok=True)
                        inj_tools.rsync(fft_file[0], cwd)
                        inj_tools.rsync(inf_file[0], cwd)
                        self.jobs.append((cwd, root, seg_args))

    def run_accelsearch(self, job):
        cwd, root, seg_args = job

        wmax = f"-wmax {seg_args['wmax']}" if seg_args.get('wmax', 0) else ''
        cmd = f"accelsearch -numharm {seg_args['numharm']} -zmax {seg_args['zmax']} {wmax} {cwd}/{root}.fft"
        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=cwd)

    def run_search(self, ncpus):
        with Pool(ncpus) as p:
            p.map(self.run_accelsearch, self.jobs)

    def transfer_products(self):
        accel_dir = f'{self.results_dir}/processing/PRESTO/ACCEL'
        os.makedirs(accel_dir, exist_ok=True)

        inj_tools.rsync(f'{self.work_dir}/*/*.inf', accel_dir)
        inj_tools.rsync(f'{self.work_dir}/*/*ACCEL_*0', accel_dir)

        fft_dir = f'{self.results_dir}/processing/PRESTO/FFT'
        if not self.processing_args['presto_search_args'].get('save_fft', False):
            for cwd, root, _ in self.jobs:
                for ext in ('.fft', '.inf'):
                    f = f'{fft_dir}/{root}{ext}'
                    if os.path.exists(f):
                        os.remove(f)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto accelsearch (stage 3/3) for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='batch process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    search_exec = PrestoAccelsearchProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    search_exec.setup()
    search_exec.run_search(args.ncpus)
    search_exec.transfer_products()
