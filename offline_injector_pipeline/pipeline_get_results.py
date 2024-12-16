import os
import json
import glob
import pandas as pd
from pathlib import Path


class ResultCollector:
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.reports = self.get_pulsar_reports()
        
    def get_pulsar_reports(self):
        reports = []
        num_dirs = sum(1 for item in os.listdir(self.output_dir)
                       if os.path.isdir(os.path.join(self.output_dir, item)) and item.startswith("inj"))

        for i in range(num_dirs):
            report_file = glob.glob(f"{self.output_dir}/inj_{i+1:06}/report_*")
            with open(report_file[0], 'r') as file:
                reports.append(json.load(file))

        return reports
    
    def get_snr(self):
        pulsars = []
        for i, report_i in enumerate(self.reports):
            directory = f'{self.output_dir}/inj_{i+1:06}/inj_pulsars'
            for pulsar in report_i['pulsars']:
                cand_file = glob.glob(f"{directory}/{pulsar['ID']}*.cands")
                cand_df = pd.read_csv(cand_file[0], skiprows=11, engine='python', delim_whitespace=True)
                PB = pulsar['binary_period'] if pulsar['binary_period'] else 0
                pulsars.append([f'inj_{i+1:06}', pulsar['ID'], pulsar['SNR'], pulsar['PX'][0], pulsar['DM'], PB, cand_df['S/N_new'].values[0]])

        collect_results = pd.DataFrame(pulsars, columns=['dir', 'ID', 'SNR', 'P0', 'DM', 'PB', 'SNR_out'])
        collect_results.to_csv(f'{self.output_dir}/0results.csv')


if __name__=='__main__':

    rc = ResultCollector('/hercules/scratch/rsenzel/offline_injections/outputs/test_1000')
    rc.get_snr()