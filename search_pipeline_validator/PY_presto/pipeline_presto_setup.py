import os
import glob
import argparse

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools


class PrestoSetup:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.inj_tag = f'inj_{self.injection_number:06}'

    def setup_processing_dir(self):
        results_dir = f'{self.out_dir}/{self.inj_tag}'
        self.processing_dir = f'{results_dir}/processing/PRESTO'
        os.makedirs(self.processing_dir, exist_ok=True)

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/{self.inj_tag}'
        report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)

    def generate_process_plan(self):
        s_args = self.processing_args['presto_search_args']

        trials = inj_tools.build_dm_trials(s_args['ddplan'], s_args.get('inj_DM', True), self.injection_report)
        batches = inj_tools.batch_trials(trials, s_args.get('batch_size', 10))

        process_tags = [f'{self.inj_tag}_BATCH_{i}' for i in range(len(batches))]

        for loc in [self.work_dir, self.processing_dir]:
            with open(f'{loc}/{self.inj_tag}_PROCESS_PLAN.txt', 'w') as f:
                for tag in process_tags:
                    f.write(tag + "\n")


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='presto setup for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    setup_exec = PrestoSetup(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    setup_exec.setup_processing_dir()
    setup_exec.get_injection_report()
    setup_exec.generate_process_plan()
