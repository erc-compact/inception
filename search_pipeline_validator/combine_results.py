import numpy as np
import pandas as pd
from pathlib import Path
import glob
import json
import os

output_dir = 'path/outputs'
inj_directories = glob.glob(f'{output_dir}/inj_*')
inj_directories.sort()
tag = 'PEASOUP'

results = []
fault = 0
keys = ['inj_dir', 'fb', 'psr_ID', 
        'inj_p0', 'inj_DM', 'inj_SNR', 'inj_DC', 'inj_s_index', 'inj_phase', 'inj_bin_p0', 'inj_M2', 'inj_inc', 'inj_ecc', 'inj_AoP', 'inj_T0',
        'par_p0', 'par_dm', 'par_acc', 'par_SNR',
        'pea_p0', 'pea_dm', 'pea_acc', 'pea_SNR', 'pea_cand_ID', 'pea_nbins', 'pea_drift', 'pea_harmonic','sift',
        'fold_p0', 'fold_dm', 'fold_acc', 'fold_SNR', 'PICS_fscore', 'PICS_recall']

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

    merged_cand_path = glob.glob(f'{inj_dir}/processing/*_SIFTED_candidates.csv')
    merged_cands = pd.read_csv(merged_cand_path[0])

    cand_paths = glob.glob(f'{inj_dir}/inj_cands/*_candfold.csv')

    pics_path = glob.glob(f'{inj_dir}/inj_cands/pics_scores.txt')
    if cand_paths:
        cand_fold = pd.read_csv(cand_paths[0], index_col=0)
        if pics_path:
            pics_vals = pd.read_csv(pics_path[0])
        else:
            pics_vals = []
    else:
        cand_fold = []
        pics_vals = []

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


        pea_results = []
        if np.any(merged_cands):
            pea_pars = merged_cands[merged_cands.get('inj_id', '') == psr['ID']]
       
            if np.any(pea_pars):
                pea_results.append(pea_pars)

        if pea_results:
            detections = pd.concat(pea_results)
            d_pea = detections.sort_values('snr', ascending=False).iloc[0]

            pea_harmonic =  int(np.round(psr['PX'][0]/d_pea['period']))
            pea_results = [d_pea['period'], d_pea['dm'], d_pea['acc'], d_pea['snr'],  
                           d_pea['cand_id'],  d_pea['nbins_offset'],  d_pea['accel_bin_drift'], pea_harmonic, d_pea['sifted']]
            inj_res.extend(pea_results)

            if len(cand_fold):
                cand_psr = cand_fold.iloc[d_pea['index_number']]
                inj_res.append(1/cand_psr['f0'])
                inj_res.append(cand_psr['dm'])
                inj_res.append(cand_psr['acc'])
                inj_res.append(cand_psr['SNR'])

                if len(pics_vals):
                    psr_pics = pics_vals.iloc[d_pea['index_number']]
                    inj_res.append(psr_pics[pics_vals.columns[1]])
                    inj_res.append(psr_pics[pics_vals.columns[2]])
                else:
                    inj_res.append(0)
                    inj_res.append(0)

            else:
                print('no cands!')
                fault += 1
                inj_res.extend(np.zeros(6).tolist())

        else:
            inj_res.extend(np.zeros(15).tolist())

        results.append(inj_res)


print(fault)
results_df = pd.DataFrame(results, columns=keys)
results_df.to_csv(f'{os.getcwd()}/results.csv')