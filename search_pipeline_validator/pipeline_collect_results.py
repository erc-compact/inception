import os
import glob
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


import pipeline_tools as inj_tools

import candidate_tools as cand_tools


def load_report(inj_dir):
    report_path = glob.glob(f'{inj_dir}/report*.json')[0]

    report = inj_tools.parse_JSON(report_path)
    pulsars = report['pulsars']
    header = report['injection_report']

    return pulsars, header


class Collector:
    def __init__(self, processing_args, out_dir, work_dir):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
    
        self.c_args = self.processing_args['collector_args']

    def load_segments(self):
        segments = []
        p_args = self.processing_args['peasoup_args']

        for s_plan in p_args['segment_plan'].keys():
            for si in range(int(s_plan)):
                segments.append(f'MATCHED_SEG_{si}_{s_plan}')

        return segments
    
    def load_pulsarx_parfold(self, report, inj_dir):
        if self.c_args['pulsarx_parfold']:
            return cand_tools.pulsarx_par2csv(report, f'{inj_dir}/inj_pulsars')
        else:
            return []
        
    def load_pulsarx_candfolds(self, inj_dir):
        segment_cands = []
        if self.c_args['pulsarx_candfold']:
            self.c_args['peasoup_search'] = True

            segments = self.load_segments()
            for segment in segments:
                cand_file = glob.glob(f'{inj_dir}/inj_cands/PEASOUP/{segment}/*.cands')
                matched_file = glob.glob(f'{inj_dir}/processing/PEASOUP/{segment}*.csv')
                if cand_file:
                    pulsarx_candfold = cand_tools.pulsarx_cand2csv(cand_file[0])
                    peasoup_df = pd.read_csv(matched_file[0], index_col=0)
                    pulsarx_candfold['segment'] = segment
                    pulsarx_candfold['PSR_ID'] = peasoup_df['PSR_ID'].values
                    segment_cands.append(pulsarx_candfold)

        if segment_cands:
            return pd.concat(segment_cands)
        else:
            return []
    
    def load_peasoup(self, inj_dir):
        segment_cands = []
        if self.c_args['peasoup_search']:
            segments = self.load_segments()
            for segment in segments:
                matched_file = glob.glob(f'{inj_dir}/processing/PEASOUP/{segment}*.csv')
                if matched_file:
                    peasoup_df = pd.read_csv(matched_file[0], index_col=0)
                    segment_cands.append(peasoup_df)

        if segment_cands:
            return pd.concat(segment_cands)
        else:
            return []
        
    def load_injection(self, psr, inj_number, max_AX, max_PX, max_FX):
        inj_results = [inj_number]
        inj_keys = ['number']

        for key, value in psr.items():

            if key == "PX":
                for i in range(max_PX):
                    inj_keys.append(f"P{i}")
                    inj_results.append(value[i] if i < len(value) else 0)

            elif key == "FX":
                for i in range(max_FX):
                    inj_keys.append(f"F{i}")
                    inj_results.append(value[i] if i < len(value) else 0)

            elif key == "AX":
                for i in range(max_AX):
                    inj_keys.append(f"A{i}")
                    inj_results.append(value[i] if i < len(value) else 0)

            elif key == "profile" and isinstance(value, dict):
                n_components = len(value["duty_cycle"])
                for i in range(n_components):
                    inj_keys.extend([f"C{i}_DC", f"C{i}_phase", f"C{i}_amp"])
                    inj_results.extend([value["duty_cycle"][i], value["phase"][i], value["amp"][i]])

            else:
                inj_keys.append(key)
                inj_results.append(value)

        inj_keys = [f'INJ_{key}' for key in inj_keys]
        return inj_keys, inj_results

    def get_input_limits(self):
        inj_directories = glob.glob(f'{self.out_dir}/inj_*')
        inj_directories.sort()

        AX, PX, FX = [], [], []
        for i, inj_dir in enumerate(inj_directories):
            report, _ = load_report(inj_dir)
            for psr in report:
                AX.append(len(psr['AX']))
                PX.append(len(psr['PX']))
                FX.append(len(psr['FX']))

        return max(AX), max(PX), max(FX), inj_directories
    
    def collect(self):
        max_AX, max_PX, max_FX, inj_directories = self.get_input_limits()

        collected_results = []
        for i, inj_dir in enumerate(inj_directories):
            inj_number = Path(inj_dir).stem
            inj_tools.print_exe(f'collecting {inj_number} ...')

            report, _ = load_report(inj_dir)

            pulsarx_parfold = self.load_pulsarx_parfold(report, inj_dir)
            pulsarx_cands = self.load_pulsarx_candfolds(inj_dir)
            peasoup_matched = self.load_peasoup(inj_dir)


            for psr in report:
                inj_keys, inj_results = self.load_injection(psr, inj_number, max_AX, max_PX, max_FX)
                
                if self.c_args['pulsarx_parfold']:
                    psr_par = pulsarx_parfold[pulsarx_parfold['PSR_ID'] == psr['ID']].iloc[0]

                    inj_results.extend([psr_par['F0'], psr_par['F0_err'], psr_par['DM'], psr_par['DM_err'], psr_par['acc'],  psr_par['acc_err'], psr_par['SNR'],  psr_par['width']])
                    parfold_keys = ['F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width']
                    inj_keys.extend([f'PAR_{key}' for key in parfold_keys])

                if self.c_args['peasoup_search']:
                    psr_pea_matched = peasoup_matched[peasoup_matched['PSR_ID'] == psr['ID']]
                    peasoup_keys = ['N_matched', 'P0', 'DM', 'acc', 'SNR', 'pepoch', 'n_samples', 'tscrunch', 'segment']

                    if len(psr_pea_matched) == 0:
                        inj_results.extend(list(np.zeros(len(peasoup_keys))))
                    else:
                        psr_pea_matched = psr_pea_matched.sort_values(by='snr', key=abs, ascending=False)
                        psr_pea = psr_pea_matched.iloc[0]
                        inj_results.extend([len(psr_pea_matched), psr_pea['period'], psr_pea['dm'], psr_pea['acc'], psr_pea['snr'], psr_pea['pepoch'], psr_pea['n_samples'], psr_pea['tscrunch'], psr_pea['segment']])
                    
                    inj_keys.extend([f'PEA_{key}' for key in peasoup_keys])

                if self.c_args['pulsarx_candfold']:
                    psr_cand_matched = pulsarx_cands[pulsarx_cands['PSR_ID'] == psr['ID']]
                    cand_keys = ['segment', 'F0', 'F0_err', 'DM', 'DM_err', 'acc', 'acc_err', 'SNR', 'width']

                    if len(psr_cand_matched) == 0:
                        inj_results.extend(list(np.zeros(len(cand_keys))))
                    else:
                        psr_cand_matched = psr_cand_matched.sort_values(by='SNR', key=abs, ascending=False)
                        psr_cand = psr_cand_matched.iloc[0]
                        inj_results.extend([psr_cand['segment'], psr_cand['F0'], psr_cand['F0_err'], psr_cand['DM'], psr_cand['DM_err'], psr_cand['acc'],  psr_cand['acc_err'], psr_cand['SNR'],  psr_cand['width']])
                    
                    inj_keys.extend([f'CAND_{key}' for key in cand_keys])

                collected_results.append(inj_results)

        results_df = pd.DataFrame(collected_results, columns=inj_keys)

        run = Path(self.out_dir).stem
        results_df.to_csv(f'{self.out_dir}/RESULTS_{run}.csv')



if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='injection results collector for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    c_exec = Collector(args.processing_args, args.out_dir, args.work_dir)
    c_exec.collect()