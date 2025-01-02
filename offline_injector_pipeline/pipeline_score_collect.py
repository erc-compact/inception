import os
import re
import glob
import json
import argparse
import subprocess
import numpy as np
import pandas as pd


from pipeline_tools import PipelineTools

class ScoreAndCollect:
    def __init__(self, pics_code, pics_models, injection_number, out_dir):
        self.pics_code = pics_code
        self.pics_models = pics_models
        self.out_dir = f'{out_dir}/inj_{injection_number:06}'
        self.inj_report = self.get_inj_report()

    def get_inj_report(self):
        report_file_path = glob.glob(f"{self.out_dir}/report_*")
        with open(report_file_path[0], 'r') as file:
            report = json.load(file)
        
        return report

    def pics_score(self, in_path):
        attempts = 0
        while (not os.path.exists(f'{in_path}/pics_scores.txt')) and (attempts != 3):
            attempts += 1
            cmd = f"python2 {self.pics_code} --in_path={in_path} --model_dir={self.pics_models} --work_dir={os.getcwd()}"
            subprocess.run(cmd, shell=True)

    def candfile_reader(self, psr_name):
        cand_file = glob.glob(f"{self.out_dir}/inj_pulsars/{psr_name}*.cands")
        cand_df = pd.read_csv(cand_file[0], skiprows=11, engine='python', sep=r'\s+')
        return cand_df

    def get_inj_results(self):
        results = pd.DataFrame(self.inj_report['pulsars'])
        results['P0'] = results['PX'].apply(lambda x: x[0])
        results = results[['ID', 'seed', 'P0', 'phase_offset', 'DM', 'SNR', 'duty_cycle', 'spectral_index',
                           'binary_period', 'T0', 'x', 'M1', 'M2', 'inc', 'ecc', 'AoP']].add_prefix('inj_')
        
        pointing_id = re.search(r'MMGPS_U_\d{4}', self.inj_report['injection']['fb']).group()
        inj_beam_name = re.search(r'[ci]fbf\d{5}', self.inj_report['injection']['fb']).group()

        results['pointing_id'] = [pointing_id] * len(results)
        results['beam_id'] = [inj_beam_name] * len(results)
        return results
    
    def get_parfold_results(self, results):
        psr_candfiles = []
        for psr in self.inj_report['pulsars']:
            cand_df = self.candfile_reader(psr['ID'])
            psr_candfiles.append(cand_df[['f0_new', 'acc_new', 'S/N_new']].add_prefix('parfold_'))
        psr_df = pd.concat(psr_candfiles)
        results[psr_df.keys()] = psr_df.values
        return results

    def get_parfold_pics_results(self, results):
        psr_par_pics = f'{self.out_dir}/inj_pulsars/pics_scores.txt'
        if os.path.exists(psr_par_pics):
            pics_par_df = pd.read_csv(psr_par_pics).drop(['arfile'], axis=1)
            pics_par_df.columns = [f'parfold_{col.split('_')[-1].split('.')[0]}' for col in pics_par_df.columns]
            results[pics_par_df.keys()] = pics_par_df.values
        return results

    def get_xmlcand_results(self, results):
        n_harm = 1
        found_files = glob.glob(f'{self.out_dir}/processing/injected_xml_candidates_harm_{n_harm}*')
        xml_cands = pd.concat([pd.read_csv(file, index_col=0) for file in found_files])
        xml_cols = ['pulsar_id', 'delta_dm', 'doppler_max', 'doppler_min', 'period', 'dm', 'acc', 'snr', '@id', 'tscrunch']
        xml_cands = xml_cands[xml_cols]

        xml_results = pd.DataFrame(columns=xml_cols)
        for i, psr in enumerate(self.inj_report['pulsars']):
            xml_psr = xml_cands[xml_cands['pulsar_id'] == psr['ID']]
            if np.any(xml_psr):
                xml_psr = xml_psr.sort_values('snr', ascending=False)
                xml_results.loc[i] = xml_psr.iloc[0]
            else:
                non_detection = np.zeros_like(xml_cols).tolist()
                non_detection[0] = psr['ID']
                xml_results.loc[i] = non_detection

        xml_results = xml_results.add_prefix('xml_')
        results[xml_results.keys()] = xml_results.values
        return results

    def get_siftcand_results(self, results):
        n_harm = 1
        sifted_csv = pd.read_csv(glob.glob(f'{self.out_dir}/inj_cands/injected_csv_candidates_harm_{n_harm}.csv')[0])
        sifted_candfile = self.candfile_reader(glob.glob(f'{self.out_dir}/inj_cands/*.cands')[0])
        sifted_candfile = sifted_candfile[:len(sifted_csv)]

        cand_cols = ['beam_index', 'tscrunch', 'adj_period']
        sifted_candfile_cols = ['dm_new', 'f0_new', 'acc_new', 'S/N_new']

        result_cols = cand_cols + sifted_candfile_cols
        cand_results = pd.DataFrame(columns=result_cols)
        for i, psr in enumerate(self.inj_report['pulsars']): 
            csv_psr = sifted_csv[sifted_csv['pulsar_id'] == psr['ID']]
            csv_candfile = sifted_candfile[sifted_csv['pulsar_id'] == psr['ID']]
            
            if np.any(csv_psr):
                csv_psr = csv_psr.sort_values('snr', ascending=False)
                csv_candfile = csv_candfile.sort_values('S/N', ascending=False)
                csv_psr = csv_psr[cand_cols]
                csv_candfile = csv_candfile[sifted_candfile_cols]
                
                combined_csv = pd.concat([csv_psr, csv_candfile], axis=1)
                cand_results.loc[i] = combined_csv.iloc[0]
            else:
                cand_results.loc[i] = np.zeros_like(result_cols)
        
        cand_results = cand_results.add_prefix('sift_')
        results[cand_results.keys()] = cand_results.values
        return results

    def get_cand_pics_results(self, results):
        psr_sift_pics = f'{self.out_dir}/inj_cands/pics_scores.txt'
        if os.path.exists(psr_sift_pics):
            n_harm = 1 
            sifted_csv = pd.read_csv(glob.glob(f'{self.out_dir}/inj_cands/injected_csv_candidates_harm_{n_harm}.csv')[0])

            pics_sift_df = pd.read_csv(psr_sift_pics).drop(['arfile'], axis=1)
            pics_cols = [f'sift_{col.split('_')[-1].split('.')[0]}' for col in pics_sift_df.columns]
            pics_sift_df.columns = pics_cols
            pics_sift_df = pics_sift_df[:len(sifted_csv)]
            
            pics_results = pd.DataFrame(columns=pics_cols)
            for i, psr in enumerate(self.inj_report['pulsars']): 
                csv_psr = sifted_csv[sifted_csv['pulsar_id'] == psr['ID']]
                pics_sift_df = pics_sift_df[sifted_csv['pulsar_id'] == psr['ID']]

                if np.any(pics_sift_df):
                    combined_df = pd.concat([pics_sift_df, csv_psr], axis=1)
                    combined_df = combined_df.sort_values('snr', ascending=False)

                    pics_results.loc[i] = combined_df[pics_cols].iloc[0]
                else:
                    pics_results.loc[i] = np.zeros_like(pics_cols)
            
            results[pics_sift_df.keys()] = pics_sift_df.values
        return results

    def collect_results(self):
        results = self.get_inj_results()
        results = self.get_parfold_results(results)
        results = self.get_parfold_pics_results(results)
        results = self.get_xmlcand_results(results)
        results = self.get_siftcand_results(results)
        results = self.get_cand_pics_results(results)

        results.to_csv(f'{self.out_dir}/inj_results.csv')
        
    def run_cmd(self):
        cand_dir = f'{self.out_dir}/inj_cands'
        par_dir = f'{self.out_dir}/inj_pulsars'

        self.pics_score(cand_dir)        
        self.pics_score(par_dir)
        
        self.collect_results()


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='PICS scorer for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--pics_code', metavar='file', required=True,  help='path to PICS code')
    parser.add_argument('--pics_models', metavar='dir', required=True, help='path to trained PICS models')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    ps = ScoreAndCollect(args.pics_code, args.pics_models, args.injection_number, args.out_dir)
    ps.run_cmd()