import re
import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Process

from pipeline_tools import PipelineTools
from pipeline_inj_cand_sifter import CandFinder
from inception.injector.io_tools import FilterbankReader



class FoldScoreExec(PipelineTools):
    def __init__(self, mode, search_args, inject_file, out_dir, injection_number, ncpus):
        super().__init__(search_args)
        self.work_dir = os.getcwd()
        self.mode = mode
        self.ncpus = ncpus
        self.out_dir = out_dir
        self.injection_number = injection_number

        filterbank, injection_report, cand_file = self.get_inputs()
        self.fb = filterbank
        self.inject_file = inject_file
        self.inj_report_path = injection_report
        self.inj_report = self.parse_JSON(injection_report)
        self.fold_cands = cand_file

        self.template = "/home/psr/software/PulsarX/include/template/meerkat_fold.template"
        self.zap_string = self.create_zap_sting()
        self.beam_tag = self.get_beam_tag()

    def get_inputs(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        for filename in os.listdir(results_dir):
            if f'_{self.data_ID}_' in filename:
                subprocess.run(f"cp {results_dir}/{filename} {self.work_dir}", shell=True)
                filterbank = f'{self.work_dir}/{Path(filename).name}'
            elif 'report' in filename:
                injection_report = f'{results_dir}/{filename}'

        if self.mode == 'cand':
            cand_file = pd.read_csv(f'{results_dir}/processing/good_cands_to_fold_with_beam.csv')
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

    def create_cand_file(self):
        self.add_tscrunch()
        self.add_adjusted_periods()

        cand_finder = CandFinder(self.fb, self.inject_file, self.inj_report_path)
        cands_df = cand_finder.parse_csv_file(self.fold_cands)
        self.n_harmonics = 2
        cands_harmon = []
        for n in range(self.n_harmonics):
            cands_data = cand_finder.filter_df(cands_df, snr_limit=5, pfact=n+1, adjust=0.1, period_key='adj_period')
            cands_data.to_csv(f'{self.work_dir}/injected_csv_candidates_harm_{n+1}.csv')
            cands_harmon.append(cands_data)

        cand_file_path = f'{self.work_dir}/candidates.candfile'
        n_cands = 0
        with open(cand_file_path, 'w') as file:
            file.write("#id DM accel F0 F1 S/N\n")
            for n in range(self.n_harmonics):
                for i, cand in cands_harmon[n].iterrows():
                    n_cands += 1
                    file.write(f"{i} {cand['dm']} {cand['acc']} {1/cand['adj_period']} 0 {cand['snr']}\n")

        return cand_file_path, n_cands
        
    def create_zap_sting(self):
        cmask = self.fold_args['channel_mask']
        cmask = cmask.strip()
        zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        return zap_string
    
    def get_beam_tag(self):
        beam_name = re.search(r'[ci]fbf\d{5}', self.inj_report['injection_report']['fb']).group()
        if 'ifbf' in beam_name:
            beam_tag = "--incoherent"
        elif 'cfbf' in beam_name:
            beam_tag = f"-i {int(beam_name.strip('cfbf'))}"
        return beam_tag

    def fold_inj_cands(self):
        cand_file_path, n_cands = self.create_cand_file()

        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)
        
        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v -t {self.ncpus} --candfile {cand_file_path} -n {nsubband} {self.beam_tag} " \
              f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {self.template} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
              f"-f {self.fb} --rfi zdot {self.zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']} -o {self.work_dir}/inj_cand"
        
        if n_cands:
            subprocess.run(cmd, shell=True)        

    def fold_par_file(self, psr):
        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)

        psr_id = psr['ID']
        par_file =  f'{self.out_dir}/inj_{self.injection_number:06}/inj_pulsars/{psr_id}.par'

        psr_p0 = psr['PX'][0]
        if psr_p0 > 0.1:
            block_size = 10
        else:
            block_size = 2

        cmd = f"psrfold_fil2 --dmboost 250 --blocksize {block_size} --plotx -v --nosearch --parfile {par_file} -n {nsubband} {self.beam_tag} " \
            f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {self.template} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
            f"-f {self.fb} --rfi zdot {self.zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']} -o {self.work_dir}/{psr_id}"

        tmp_cwd = f'{self.work_dir}/process_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)
        subprocess.run(cmd, shell=True, cwd=tmp_cwd) 

    def fold_par_pulsars(self):

        processes = []
        for psr in self.inj_report['pulsars']:
            process = Process(target=self.fold_par_file, args=(psr, ))
            processes.append(process)
            process.start()
        
        for process in processes:
            process.join()    
    
    def delete_inj_filterbank(self):
        fb_name = Path(self.fb).name
        fb_path = f'{self.out_dir}/inj_{self.injection_number:06}/{fb_name}'
        os.remove(fb_path)

    def run_cmd(self):
        if self.mode == 'par':
            self.fold_par_pulsars()
        elif self.mode == 'cand':
            self.fold_inj_cands()
            self.delete_inj_filterbank()

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        if self.mode == 'par':
            par_dir = f'{results_dir}/inj_pulsars'
            subprocess.run(f"cp {self.work_dir}/*.png {par_dir}", shell=True)
            subprocess.run(f"cp {self.work_dir}/*.ar {par_dir}", shell=True)
            subprocess.run(f"cp {self.work_dir}/*.cands {par_dir}", shell=True)

        elif self.mode == 'cand':
            cand_dir = f'{results_dir}/inj_cands'
            if os.path.isfile(cand_dir):
                os.remove(cand_dir)
            os.makedirs(cand_dir, exist_ok=True)
            subprocess.run(f"cp {self.work_dir}/*.png {cand_dir}", shell=True)
            subprocess.run(f"cp {self.work_dir}/*.ar {cand_dir}", shell=True)
            subprocess.run(f"cp {self.work_dir}/*.cands {cand_dir}", shell=True)
            subprocess.run(f"cp {self.work_dir}/candidates.candfile {cand_dir}", shell=True)
            for n in range(self.n_harmonics):
                subprocess.run(f"cp {self.work_dir}/injected_csv_candidates_harm_{n+1}.csv {cand_dir}", shell=True)
                

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate folder for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--mode', metavar='str', required=True,  help='which mode to fold')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_file', metavar='file', required=True, help='JSON file with injection plan')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    parser.add_argument('--ncpus', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = FoldScoreExec(args.mode, args.search_args, args.injection_file, args.out_dir, args.injection_number, args.ncpus)
    fold_exec.run_cmd()
    fold_exec.transfer_products()
   