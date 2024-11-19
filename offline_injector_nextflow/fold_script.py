import re
import sys
import json
import argparse
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Pool
from collections import namedtuple

from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

from inception.injector.io_tools import FilterbankReader


class FoldScoreExec:
    def __init__(self, fb, search_args, injection_report, fold_cands, output, par_dir, num_threads):
        self.out = output
        self.par_dir = par_dir
        self.fb = fb
        self.num_threads = num_threads

        args = self.parse_JSON(search_args)
        self.fold_args = args['fold_args']
        self.data_args = args['data']
        self.processing_args = args['processing_args']
        self.inj_report = self.parse_JSON(injection_report)
        self.fold_cands = pd.read_csv(fold_cands)
    
    @staticmethod
    def parse_fb(fb):
        if type(fb) == list:
            return ' '.join(fb)
        else:
            return fb

    def parse_JSON(self, json_file):
        try:
            with open(json_file, 'r') as file:
                pars = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {json_file}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {json_file} using JSON.')
        else:
            return pars
        
    def create_DDplan(self):
        DMRange = namedtuple("DMRange", ["low_dm", "high_dm", "dm_step", "tscrunch"])

        segments = []
        plan = self.processing_args['ddplan']
        for line in plan.splitlines():
            low_dm, high_dm, dm_step, tscrunch = list(map(float, line.split()[:4]))
            segments.append(DMRange(low_dm, high_dm, dm_step, tscrunch))

        return list(sorted(segments, key=lambda x: x.tscrunch))
            
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

        cand_file_path = f'{self.out}/candidates.candfile'
        with open(cand_file_path, 'w') as file:
            file.write("#id DM accel F0 F1 S/N\n")
            for i, cand in cands_data.iterrows():
                file.write(f"{i} {cand['dm']} {cand['acc']} {1/cand['adj_period']} 0 {cand['snr']}\n")
        file.close()

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
        template="/home/psr/software/PulsarX/include/template/meerkat_fold.template"
        zap_string = self.create_zap_sting()
        beam_tag = self.get_beam_tag()
        cand_file = self.create_cand_file()

        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)
        
        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v -t {self.num_threads} --candfile {cand_file} -n {nsubband} {beam_tag} " \
              f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {template} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
              f"-f {self.fb} --rfi zdot {zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']} -o {self.out}/inj_cand"
    
        subprocess.run(cmd, shell=True)

    def pics_score(self):
        # pics_code = '/home/psr/software/trapum-pipeline-wrapper/pipelines/trapum_fold_and_score_pipeline/webpage_score.py'
        pics_code = '/u/rsenzel/pulsar_inject/nextflow/webpage_score_hercules.py'
        cmd = f"python2 {pics_code} --in_path={self.out}"
        subprocess.run(cmd, shell=True)

    def fold_pars(self):

        def fold(psr_id):
            if psr_id:
                par_file = f"{self.par_dir}/{psr_id}.par"

                cmd = f"psrfold_fil2 --dmboost 250 --plotx -v --parfile {par_file} -n {nsubband} {beam_tag} " \
                    f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {template} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
                    f"-f {self.fb} --rfi zdot {zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']} -o {self.out}/{psr_id}"

                subprocess.run(cmd, shell=True)  

        template="/home/psr/software/PulsarX/include/template/meerkat_fold.template"
        zap_string = self.create_zap_sting()
        beam_tag = self.get_beam_tag()

        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)

        psr_ids = [psr_id['ID'] for psr_id in self.inj_report['pulsars']]

        args_list = [[] for _ in range(self.num_threads)]
        for i, file in enumerate(psr_ids):
            args_list[i % self.num_threads].append(file)

        with Pool(self.num_threads) as p:
            p.map(fold, args_list)

    def collect_results(self):
        pass


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate folder for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--fb', metavar='file', required=True,  help='filterbank(s) to fold')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_report', metavar='file', required=True, help='JSON file with inject pulsar records')
    parser.add_argument('--fold_cands', metavar='file', required=True, help='csv file with good_cands_to_fold_with_beam')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    parser.add_argument('--par_dir', metavar='dir', required=True, help='directory with injected pulsar parfiles')
    parser.add_argument('--n_threads', metavar='int', type=int, required=True, help='number of threads to use')
    args = parser.parse_args()

    fold_exec = FoldScoreExec(args.fb, args.search_args, args.injection_report, args.fold_cands, args.output, args.par_dir, args.n_threads)
    fold_exec.fold_inj_cands()
    fold_exec.pics_score()
    fold_exec.fold_pars()



