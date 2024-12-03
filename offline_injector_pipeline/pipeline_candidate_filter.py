import os
import re
import argparse
import subprocess
import pandas as pd
from pathlib import Path

from pipeline_tools import PipelineTools


class CandExec(PipelineTools):
    def __init__(self, search_args, data_dir, out_dir, injection_number):
        super().__init__(search_args)
        self.work_dir = os.getcwd()
        self.data_dir = data_dir
        self.out_dir = out_dir
        self.injection_number = injection_number

        self.inj_report = self.get_injection_report()
        self.injection_ID = self.inj_report['injection']['ID']
    
    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        for filename in os.listdir(results_dir):
            if 'report' in filename:
                return self.parse_JSON(filename)
    
    def create_XML_list(self):
        ddplan = self.create_DDplan()
        xml_file_names = [f'overview_dm_{dm_range.low_dm:.6f}_{dm_range.high_dm:.6f}.xml' for dm_range in ddplan]

        pointing_id = re.search(r'MMGPS_U_\d{4}', self.inj_report['injection']['fb']).group()
        inj_beam_name = re.search(r'[ci]fbf\d{5}', self.inj_report['injection']['fb']).group()
        xml_file_paths = []
        for nbeam in range(self.data['n_cbeams']):
            for xml_name in xml_file_names:
                beam_i = f'cfbf{nbeam:05g}'
                if beam_i == inj_beam_name:
                    xml_file_paths.append(f'{self.out_dir}/inj_{self.injection_number:06}/processing/{self.data_ID}_{self.injection_ID}_{xml_name}')
                else:
                    xml_file_paths.append(f'{self.data_dir}/{pointing_id}/XML_FILES/{beam_i}/{xml_name}')
        if self.data['n_ibeams']:
            for xml_name in xml_file_names:
                xml_file_paths.append(f'{self.data_dir}/{pointing_id}/XML_FILES/ifbf00000/{xml_name}')

        outfile = f'{self.work_dir}/candidates_{self.data_ID}_{self.injection_ID}.ascii'
        with open(outfile, 'w') as file:
            file.writelines(line + '\n' for line in xml_file_paths)
        
        return outfile

    def run_cmd(self):
        candidate_files = self.create_XML_list()
        
        cmd = f"candidate_filter.py -i {candidate_files} -o {self.work_dir}/ --threshold {self.multi_beam['snr_cutoff']} " \
              f"--p_tol {self.multi_beam['p_tol']} --dm_tol {self.multi_beam['dm_tol']} " \
              "-c /home/psr/software/candidate_filter/candidate_filter/default_config.json " \
              "--rfi /home/psr/software/candidate_filter/candidate_filter/known_rfi.txt"

        subprocess.run(cmd, shell=True)

    def transfer_products(self):
        def get_beam_name(file_path):
            path_obj = Path(file_path)
            if path_obj.stem.split('_')[0] == 'data':
                return re.search(r'[ci]fbf\d{5}', self.inj_report['injection']['fb']).group()
            else:
                return path_obj.parent.stem 

        fold_candidates = pd.read_csv(f'{self.work_dir}/_good_cands_to_fold.csv')
        fold_candidates['beam_id'] = [get_beam_name(file_path) for file_path in fold_candidates['file']]
        
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}/processing'
        fold_candidates.to_csv(f'{results_dir}/good_cands_to_fold_with_beam.csv')
    

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate filter for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--data_dir', metavar='dir', required=True, help='directory of xml files')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    cand_exec = CandExec(args.search_args, args.data_dir, args.out_dir, args.injection_number)
    cand_exec.run_cmd()
    cand_exec.transfer_products()