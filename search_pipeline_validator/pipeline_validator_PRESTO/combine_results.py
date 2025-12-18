import numpy as np
import pandas as pd
from pathlib import Path
import glob
import json
import os

output_dir = 'path/outputs'
inj_directories = glob.glob(f'{output_dir}/inj_*')
inj_directories.sort()

results = []
fault = 0
keys = ['inj_dir', 'fb', 'psr_ID', 
        'inj_p0', 'inj_DM', 'inj_SNR', 'inj_DC', 'inj_s_index', 'inj_phase', 'inj_bin_p0', 'inj_M2', 'inj_inc', 'inj_ecc', 'inj_AoP', 'inj_T0',
        'par_p0', 'par_dm', 'par_acc', 'par_SNR', 'found']

for i, inj_dir in enumerate(inj_directories):

    inj_number = Path(inj_dir).stem
    print(inj_number)

    report_path = glob.glob(f'{inj_dir}/report*.json')[0]
    with open(report_path, 'r') as file:
        load_report = json.load(file)
        report = load_report['pulsars']
        fb_name = Path(load_report['injection_report']['fb']).stem

    par_fold_path = glob.glob(f'{inj_dir}/inj_pulsars/*_parfold.csv')[0]
    par_fold = pd.read_csv(par_fold_path, index_col=0)

    presto_df = pd.read_csv(f'{inj_dir}/inj_cands_PRESTO/candidates.csv')

    for psr in report:
        inj_res = [inj_number, fb_name, psr['ID'], psr['PX'][0], psr['DM'], psr['SNR'], psr['duty_cycle'], psr['spectral_index'], psr['phase_offset']]

        if psr['binary_period']:
            inj_res.extend([psr['binary_period'], psr['M2'], psr['inc'], psr['ecc'], psr['AoP'], psr['T0']])
        else:
            inj_res.extend(np.zeros(6).tolist())
                   
        psr_par = par_fold[par_fold['ID'] == psr['ID']].iloc[0]
        inj_res.append(1/psr_par['f0'])
        inj_res.append(psr_par['dm'])
        inj_res.append(psr_par['acc'])
        inj_res.append(psr_par['SNR'])


        df_f = presto_df[presto_df['ID'] == psr['ID']]
        if df_f:
            inj_res.extend(['found'])
        else:
            inj_res.extend(['not found'])

        results.append(inj_res)


print(fault)
results_df = pd.DataFrame(results, columns=keys)
results_df.to_csv(f'{os.getcwd()}/PRESTO_results.csv')