import os
import argparse
import subprocess

class PicsScorer:
    def __init__(self, pics_code, pics_models, injection_number, out_dir):
        self.pics_code = pics_code
        self.pics_models = pics_models
        self.injection_number = injection_number
        self.out_dir = out_dir

    def pics_score(self, in_path):
        cmd = f"python2 {self.pics_code} --in_path={in_path} --model_dir={self.pics_models} --work_dir={os.getcwd()}"
        subprocess.run(cmd, shell=True)        

    def run_cmd(self):

        par_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars'
        if os.path.isdir(par_dir):
            self.pics_score(par_dir)
        
        cand_dir = f'{self.out_dir}/inj_{self.injection_number:06}/inj_cands'
        if os.path.isdir(cand_dir):
            self.pics_score(cand_dir)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='PICS scorer for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--pics_code', metavar='file', required=True,  help='path to PICS code')
    parser.add_argument('--pics_models', metavar='dir', required=True, help='path to trained PICS models')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    ps = PicsScorer(args.pics_code, args.pics_models, args.injection_number, args.out_dir)
    ps.run_cmd()