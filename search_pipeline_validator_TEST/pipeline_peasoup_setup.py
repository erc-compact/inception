import os
import argparse

import search_pipeline_validator_TEST.pipeline_tools as inj_tools


class PeasoupSetup:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number
        self.inj_tag = f'inj_{self.injection_number:06}'

    def setup_processing_dir(self):
        results_dir = f'{self.out_dir}/{self.inj_tag}'
        self.processing_dir = f'{results_dir}/processing/PEASOUP'
        os.makedirs(self.processing_dir, exist_ok=True)

    def generate_process_plan(self):
        p_args = self.processing_args['peasoup_args']

        process_tags = []
        for d_plan in p_args['ddplan'].keys():
            for s_plan in p_args['segment_plan'].keys():
                for si in range(int(s_plan)):
                    process_tags.append(f'{self.inj_tag}_DDPLAN_{d_plan}_SEG_{si}_{s_plan}')

        for loc in [self.work_dir, self.processing_dir]:
            with open(f'{loc}/{self.inj_tag}_PROCESS_PLAN.txt', 'w') as f:
                for tag in process_tags:
                    f.write(tag + "\n")


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='peasoup setup for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    peasoup_exec = PeasoupSetup(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    peasoup_exec.setup_processing_dir()
    peasoup_exec.generate_process_plan()