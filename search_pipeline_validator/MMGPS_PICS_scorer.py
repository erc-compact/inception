import os
import glob
import argparse
import subprocess


import pipeline_tools as inj_tools


class PICSScorer:
    def __init__(self,  processing_args, injection_number, out_dir, work_dir):
        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.pics_models = self.processing_args['PICS_scorer']['models_dir']
        self.pics_code = self.processing_args['PICS_scorer']['code_path']
        self.injection_number = injection_number

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

    def scorer_setup(self):
        self.get_injection_report()
    
    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        self.report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(self.report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def run_cmd(self):
        attempts = 0
        in_path = f'{self.out_dir}/inj_{self.injection_number:06}/inj_cands'
        ar_files = glob.glob(f'{in_path}/*.ar')
        if len(ar_files) != 0:
            while (not os.path.exists(f'{in_path}/pics_scores.txt')) and (attempts != 3):
                attempts += 1
                cmd = f"python2 {self.pics_code} --in_path={in_path} --model_dir={self.pics_models} --work_dir={os.getcwd()}"
                subprocess.run(cmd, shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='MMGPS PICS scorer for INCEPTION',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--processing_args', metavar='dir', required=True, help='path to trained PICS models')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')
    args = parser.parse_args()

    ps = PICSScorer(args.processing_args, args.injection_number, args.out_dir,  args.work_dir)
    ps.scorer_setup()
    ps.run_cmd()