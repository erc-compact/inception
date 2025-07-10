import os
import glob
import argparse
import subprocess
import pandas as pd
from pathlib import Path

import pipeline_tools as inj_tools


class CandidateFilterProcess:
    def __init__(self,  processing_args, out_dir, work_dir, injection_number):
        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def filter_setup(self):
        self.get_injection_report()
        self.get_observation()
        self.create_XML_list()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        self.report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(self.report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def get_observation(self):
        fb_inj_path = self.injection_report['injection_report']['fb']
        obs_id = '_'.join(Path(fb_inj_path).stem.split('_')[:2])
        self.mjd, self.beam = obs_id.split('_')
        data = self.processing_args['data']

        beam_path = []
        for fb in data:
            if type(fb) == list:
                fb_path = fb[0]
            else:
                fb_path = fb
            if obs_id in fb_path:
                beam_path.append(fb_path)

        self.xml_dir = f'{str(Path(beam_path[0]).parent.parent)}/XML_FILES'

    def create_XML_list(self):
        ddplan = inj_tools.create_DDplan(self.processing_args['peasoup_args']['ddplan'])
        xml_file_names = [f'{dm_range.low_dm:.6f}_{dm_range.high_dm:.6f}.xml' for dm_range in ddplan]
        file_id = f"{self.processing_args['injection_args']['id']}_{self.inj_id}"

        n_cbeams = self.processing_args['MMGPS_candidate_filter']['n_cbeams']
        n_ibeams = self.processing_args['MMGPS_candidate_filter']['n_ibeams']
        xml_file_paths = []
        for nbeam in range(n_cbeams):
            for xml_name in xml_file_names:
                beam_i = f'cfbf{nbeam:05g}'
                if beam_i == self.beam:
                    xml_file = f"{self.out_dir}/inj_{self.injection_number:06}/processing/{file_id}_DM_{xml_name}"
                    if glob.glob(xml_file):
                        xml_file_paths.append(xml_file)
                    else:
                        xml_file_paths.append(f'{self.xml_dir}/{beam_i}/overview_dm_{xml_name}')
                else:
                    xml_file_paths.append(f'{self.xml_dir}/{beam_i}/overview_dm_{xml_name}')
        
        if n_ibeams:
            for xml_name in xml_file_names:
                xml_file_paths.append(f'{self.xml_dir}/ifbf00000/overview_dm_{xml_name}')

        self.xml_file = f'{self.work_dir}/candidates_{file_id}.ascii'
        with open(self.xml_file, 'w') as file:
            file.writelines(line + '\n' for line in xml_file_paths)

    def run_cmd(self):
        params = self.processing_args['MMGPS_candidate_filter']
        cmd = f"candidate_filter.py -i {self.xml_file} -o {self.work_dir}/ --threshold {params['snr_cutoff']} " \
              f"--p_tol {params['p_tol']} --dm_tol {params['dm_tol']} " \
              "-c /home/psr/software/candidate_filter/candidate_filter/default_config.json " \
              "--rfi /home/psr/software/candidate_filter/candidate_filter/known_rfi.txt"

        subprocess.run(cmd, shell=True)

    def transfer_products(self):
        def get_beam_name(file_path):
            path_obj = Path(file_path)
            if path_obj.stem.split('_')[0] != 'overview':
                return self.beam
            else:
                return path_obj.parent.stem 

        fold_candidates = pd.read_csv(f'{self.work_dir}/_good_cands_to_fold.csv')
        fold_candidates['beam_id'] = [get_beam_name(file_path) for file_path in fold_candidates['file']]
        
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/processing'
        fold_candidates.to_csv(f'{results_dir}/good_cands_to_fold_with_beam.csv')



if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='MMGPS candidate filter for INCEPTION',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')
    args = parser.parse_args()

    cand_exec = CandidateFilterProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    cand_exec.filter_setup()
    cand_exec.run_cmd()
    cand_exec.transfer_products()