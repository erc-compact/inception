import re
import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Process

from pipeline_tools import PipelineTools
from inception.injector.io_tools import FilterbankReader


class FoldScoreExec(PipelineTools):
    def __init__(self, mode, search_args, out_dir, injection_number, ncpus):
        super().__init__(search_args)
        self.work_dir = os.getcwd()
        self.mode = mode
        self.ncpus = ncpus
        self.out_dir = out_dir
        self.injection_number = injection_number

        filterbank, injection_report, cand_file = self.get_inputs()
        self.fb = filterbank
        self.inj_report = self.parse_JSON(injection_report)
        self.fold_cands = cand_file

        self.template = "/home/psr/software/PulsarX/include/template/meerkat_fold.template"
        self.zap_string = self.create_zap_sting()
        self.beam_tag = self.get_beam_tag()

    def get_inputs(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        for filename in os.listdir(results_dir):
            if f'_{self.data_ID}_' in filename:
                subprocess.run(f"rsync -Pav {filename} {self.work_dir}", shell=True)
                filterbank = f'{self.work_dir}/{Path(filename).name}'
            elif 'report' in filename:
                injection_report = filename

        if self.mode == 'cand':
            cand_file = pd.read_csv(f'{results_dir}/good_cands_to_fold_with_beam.csv')
        else:
            cand_file = ''
        
        return filterbank, injection_report, cand_file
            
    def cand_cutoffs(self):

        def period_parse_cuts(cuts, tobs):
            if ":" not in cuts:
                return float(cuts)
            else:
                for cut in cuts.split(","):
                    low, high, value = list(map(float, cut.split(":")))
                    if tobs >= low and tobs < high:
                        return value
                else:
                    return 0.0
                
        fbreader = FilterbankReader(self.fb)
        period, dm, snr = self.fold_cands['period'], self.fold_cands['dm'], self.fold_cands['snr']

        tobs = fbreader.n_samples * fbreader.dt
        period_cutoff = self.fold_args.get('period_cutoffs', '0:inf:0.0000000001')
        period_cut = period_parse_cuts(period_cutoff, tobs)

        ftop, fbottom, nchans = fbreader.ftop, fbreader.fbottom, fbreader.nchans
        smearing_cutoff = np.abs(4.148741601e3 * dm * (1/(fbottom**2) - 1/(ftop**2))/nchans)

        snr_cutoff = float(self.fold_args['snr_cutoff'])
        del fbreader
        return (period > period_cut) & (period > smearing_cutoff) & (snr > snr_cutoff)
    
    def add_tscrunch(self):
        ddplan = self.create_DDplan()
        dm2tsrunch = {f'{dm_range.low_dm:.6f}':int(dm_range.tscrunch) for dm_range in ddplan}
        xml_start_dm = [Path(file).stem.split('_')[-2] for file in self.fold_cands['file']]

        self.fold_cands['tscrunch'] = [dm2tsrunch[dm_i] for dm_i in xml_start_dm]
    
    def add_adjusted_periods(self):

        def period_obs_centre(p0, pdot, tsamp, n_samples, fft_size):
            return p0 - pdot * (fft_size - n_samples) * tsamp / 2
            
        period, acc, tscrunch = self.fold_cands['period'], self.fold_cands['acc'], self.fold_cands['tscrunch']
        pdot = period * acc / 2.99792458e8

        fbreader = FilterbankReader(self.fb)
        tsamp = float(fbreader.dt)
        n_samples = int(fbreader.n_samples)
        fftsize = int(self.processing_args['fft_length'])
        self.fold_cands['adj_period'] = period_obs_centre(period, pdot, float(tsamp), n_samples//tscrunch, fftsize//tscrunch)
        del fbreader

    def filter_cands(self):
        cands_cut = self.fold_cands[self.cand_cutoffs()]
        cands_list = cands_cut.sort_values(by='snr', ascending=False)
        cands_list.index = np.arange(len(cands_list))

        max_cands = self.fold_args['cand_limit_per_beam']
        beams = cands_list['beam_id'].unique()
        
        keep_cands_index = []
        beam_counts = np.zeros_like(beams)
        for i, cand in cands_list.iterrows():
            beam_index = np.where(beams == cand['beam_id'])[0]
            if beam_counts[beam_index] != max_cands:
                keep_cands_index.append(i)
                beam_counts[beam_index] += 1

        filtered_cands = cands_list.iloc[keep_cands_index]
        filtered_cands.index = np.arange(len(filtered_cands))
        return filtered_cands
    
    def get_injected_cands(self, filtered_cands):
        psr = self.inj_report['pulsars']
        ptol, dmtol = 0.01, 0.01

        pulsar_cands = [[] for _ in range(len(psr))]
        for i, psr_i in enumerate(psr):
            p0_i = psr_i['PX'][0]
            dm_i = psr_i['DM']
            p_cond = (filtered_cands['period'] > p0_i*(1-ptol)) & (filtered_cands['period'] < p0_i*(1+ptol))
            dm_cond = (filtered_cands['dm'] > dm_i*(1-dmtol)) & (filtered_cands['dm'] < dm_i*(1+dmtol))
            psr_cands = filtered_cands[p_cond & dm_cond]
            pulsar_cands[i].extend(psr_cands.index.values)

        return filtered_cands.iloc[np.concatenate(pulsar_cands).astype(int)]

    def create_cand_file(self):
        self.add_tscrunch()
        self.add_adjusted_periods()
        filtered_cands = self.filter_cands()

        cands_data = self.get_injected_cands(filtered_cands)

        cand_file_path = f'{self.work_dir}/candidates.candfile'
        with open(cand_file_path, 'w') as file:
            file.write("#id DM accel F0 F1 S/N\n")
            for i, cand in cands_data.iterrows():
                file.write(f"{i} {cand['dm']} {cand['acc']} {1/cand['adj_period']} 0 {cand['snr']}\n")

        return cand_file_path
        
    def create_zap_sting(self):
        cmask = self.fold_args['channel_mask']
        cmask = cmask.strip()
        zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        return zap_string
    
    def get_beam_tag(self):
        beam_name = re.search(r'[ci]fbf\d{5}', self.inj_report['injection']['fb']).group()
        if 'ifbf' in beam_name:
            beam_tag = "--incoherent"
        elif 'cfbf' in beam_name:
            beam_tag = f"-i {int(beam_name.strip('cfbf'))}"
        return beam_tag

    def fold_inj_cands(self):
        cand_file = self.create_cand_file()

        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)
        
        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v -t {self.ncpus} --candfile {cand_file} -n {nsubband} {self.beam_tag} " \
              f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {self.template} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
              f"-f {self.fb} --rfi zdot {self.zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']} -o {self.work_dir}/inj_cand"

        subprocess.run(cmd, shell=True)

    def fold_par_file(self, psr_id):
        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)

        par_file =  f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr_id}.par'

        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v --nosearch --parfile {par_file} -n {nsubband} {self.beam_tag} " \
            f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {self.template} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
            f"-f {self.fb} --rfi zdot {self.zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']} -o {self.work_dir}/{psr_id}"

        tmp_cwd = f'{self.work_dir}/process_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)
        subprocess.run(cmd, shell=True, cwd=tmp_cwd) 

    def fold_par_pulsars(self):
        fold_list = [psr_id['ID'] for psr_id in self.inj_report['pulsars']]

        processes = []
        for psr_id in fold_list:
            process = Process(target=self.fold_par_file, args=(psr_id, ))
            processes.append(process)
            process.start()
        
        for process in processes:
            process.join()

    def pics_score(self):
        # pics_code = '/home/psr/software/trapum-pipeline-wrapper/pipelines/trapum_fold_and_score_pipeline/webpage_score.py'
        pics_code = '/hercules/scratch/rsenzel/offline_injections/injections_config/test_1000/PICS_MMGPS.py'
        model_dir = '/hercules/scratch/rsenzel/offline_injections/injections_config/test_1000/pics_models'
        cmd = f"python2 {pics_code} --in_path={self.work_dir} --model_dir={model_dir}"
        subprocess.run(cmd, shell=True)        

    def run_cmd(self):
        if self.mode == 'par':
            self.fold_par_pulsars()
        elif self.mode == 'cand':
            self.fold_inj_cands()

        self.pics_score()

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        if self.mode == 'par':
            par_dir = f'{results_dir}/inj_pulsars'
            subprocess.run(f"rsync -Pav {self.work_dir}/*.png {par_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/*.ar {par_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/*.cands {par_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/pics_scores.txt {par_dir}", shell=True)

        elif self.mode == 'cand':
            cand_dir = f'{results_dir}/inj_cands'
            os.mkdir(cand_dir)
            subprocess.run(f"rsync -Pav {self.work_dir}/*.png {cand_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/*.ar {cand_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/*.cands {cand_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/*.candfile {cand_dir}", shell=True)
            subprocess.run(f"rsync -Pav {self.work_dir}/pics_scores.txt {cand_dir}", shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate folder for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--mode', metavar='str', required=True,  help='which mode to fold')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    parser.add_argument('--ncpus', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = FoldScoreExec(args.mode, args.search_args, args.out_dir, args.injection_number, args.ncpus)
    fold_exec.run_cmd()
    fold_exec.transfer_products()
