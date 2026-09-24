import glob
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

import pipeline_tools as inj_tools
import candidate_tools as cand_tools

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from injector.io_tools import FilterbankReader, print_exe
from injector.setup_manager import SetupManager



class CandidateMatcher:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'


        self.candidate_loaders = {
            'PEASOUP': self.load_peasoup_candidates,
            'PRESTO': self.load_presto_candidates,
        }
        self.fold_file_generators = {
            'PEASOUP': self.generate_peasoup_fold_files,
            'PRESTO': self.generate_presto_fold_files,
        }

    def setup(self):
        self.get_injection_report()

        fb_path = glob.glob(f"{self.results_dir}/*_{self.inj_id}.fil")[0]

        ephem = self.processing_args['injection_args']['ephem']
        if ephem != 'builtin':
            inj_tools.rsync(ephem, self.work_dir)
            ephem = f'./{Path(ephem).name}'

        self.fb = FilterbankReader(fb_path, stats_samples=0)
        self.setup_manager = SetupManager(self.report_path, fb_path, ephem, generate=False, override_length=0)

    def get_injection_report(self):
        self.report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(self.report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def get_mode(self):
        return self.processing_args['candidate_matcher_args']['mode']


    @staticmethod
    def parse_peasoup_tag(xml_tag):
        splits = xml_tag.split('_')
        return splits[-2], splits[-1], splits[-4]

    def correct_freq(self, csv_cands, fftsize, dt):
        sys.exit(1)
        # cand_tools.correct_fftsize_offset(period, acc, fftsize, nsamples, dt)

    def load_peasoup_candidates(self):
        processing_dir = f'{self.results_dir}/processing'
        ref_dt = self.setup_manager.pulsar_models[0].obs.dt

        xmls = glob.glob(f'{processing_dir}/PEASOUP/*.xml')
        cands_list = []
        for xml in xmls:
            csv_cands, pepoch, fftsize, dt = cand_tools.xml2csv(xml)
            if pepoch is not None:
                csv_cands['pepoch'] = pepoch
            else:
                self.correct_freq(csv_cands, fftsize, dt)
            csv_cands['n_samples'] = fftsize
            csv_cands['T'] = fftsize * ref_dt
            csv_cands['XML'] = Path(xml).name
            s0, s1, tscrunch = self.parse_peasoup_tag(Path(xml).stem)
            csv_cands['tscrunch'] = tscrunch
            csv_cands['segment'] = f'{s0}_{s1}'
            cands_list.append(csv_cands)

        candidates = pd.concat(cands_list, ignore_index=True) if cands_list else pd.DataFrame()

        if len(candidates):
            candidates.to_csv(f'{processing_dir}/PEASOUP/inj_{self.injection_number:06}_PEASOUP_candidates.csv')
        return candidates

    def generate_peasoup_fold_files(self):
        processing_dir = f'{self.results_dir}/processing/PEASOUP'
        unique_segments = np.unique(self.comb_matched_cands['segment'])

        process_tags = []
        for segment in unique_segments:
            s0, s1 = segment.split('_')
            segment_cands = self.comb_matched_cands[self.comb_matched_cands['segment'] == segment]
            pepoch = segment_cands['pepoch'].values[0]

            process_tag = f'MATCHED_SEG_{s0}_{s1}'
            segment_cands.to_csv(f'{processing_dir}/{process_tag}_{pepoch}.csv')

            candfile_path = f'{processing_dir}/{process_tag}_{pepoch}.candfile'
            cand_tools.create_PULSARX_candfile(segment_cands, candfile_path)

            process_tags.append(process_tag)

        self.write_fold_plan(process_tags, processing_dir)

    def load_presto_candidates(self):
        processing_dir = f'{self.results_dir}/processing/PRESTO'
        sift_csv = glob.glob(f'{processing_dir}/PRESTO_candidates.csv')
        if not sift_csv:
            print_exe('No sifted PRESTO candidates found.')
            return pd.DataFrame()

        candidates = cand_tools.presto_sift2csv(sift_csv[0])

        ref_obs = self.setup_manager.pulsar_models[0].obs
        mid_time_sec = (candidates['seg_i'] + 0.5) / candidates['seg_n'] * ref_obs.obs_len
        candidates['pepoch'] = ref_obs.sec2mjd(mid_time_sec.values)

        candidates.to_csv(f'{processing_dir}/inj_{self.injection_number:06}_PRESTO_candidates.csv')
        return candidates

    def generate_presto_fold_files(self):
        processing_dir = f'{self.results_dir}/processing/PRESTO'
        unique_segments = np.unique(self.comb_matched_cands['segment'])

        process_tags = []
        for segment in unique_segments:
            segment_cands = self.comb_matched_cands[self.comb_matched_cands['segment'] == segment]

            process_tag = f'MATCHED_SEG_{segment}'
            segment_cands.to_csv(f'{processing_dir}/{process_tag}.csv')

            process_tags.append(process_tag)

        self.write_fold_plan(process_tags, processing_dir)

    def write_fold_plan(self, process_tags, processing_dir):
        for loc in [self.work_dir, processing_dir]:
            with open(f'{loc}/inj_{self.injection_number:06}_FOLD_PLAN.txt', 'w') as f:
                for tag in process_tags:
                    f.write(tag + "\n")

    def process(self):
        mode = self.get_mode()
        harmonics = self.processing_args['candidate_matcher_args'].get('harmonics', [0.5, 1, 2])

        loader = self.candidate_loaders.get(mode)
        if loader is None:
            sys.exit(f"Unknown candidate_matcher_args mode '{mode}'.")

        candidates = loader()
        if len(candidates) == 0:
            self.matched_cands = []
            self.comb_matched_cands = pd.DataFrame()
            return

        self.matched_cands = self.match_candidates(candidates, harmonics)
        self.comb_matched_cands = pd.concat(self.matched_cands, ignore_index=True) if self.matched_cands else pd.DataFrame()

    def get_obs_params(self):
        ref_psr = self.setup_manager.pulsar_models[0]
        obs = ref_psr.obs

        topo_sec = np.linspace(0, obs.obs_len, 4800)
        topo_mjd = obs.sec2mjd(topo_sec)
        bary_sec = obs.topo2bary(topo_mjd, return_mjd=False, interp=False)

        rv = obs.earth_radial_velocity(topo_mjd)
        return rv, topo_sec, bary_sec

    def match_candidates(self, candidates, harmonics):
        matched_pulsar_cands = []
        max_folds = self.processing_args['candidate_matcher_args']['max_folds']
        obs_rv, topo_sec, bary_sec = self.get_obs_params()
        bin_tol = 2

        for pm in self.setup_manager.pulsar_models:
            print_exe(f'Matching PSR {pm.ID} ...')
            fft_bin = 1 / candidates['T']
            F0_cands = 1/candidates['period']

            frame = pm.pulsar_pars['frame']
            if frame == 'bary':
                time = bary_sec
            elif frame == 'topo':
                time = topo_sec
                obs_rv = np.zeros_like(obs_rv)

            rv_PSR = cand_tools.add_PSR_rv_curve(pm, time, obs_rv.copy(), pm.pulsar_pars['ACCEPOCH']) # acc dependent, low freq psr 0 accel test

            F_min, F_max = cand_tools.get_freq_bounds(rv_PSR, pm, candidates['pepoch'])

            freq_cond = np.zeros(len(candidates), dtype=bool)
            for harmonic in harmonics:
                freq_cond |= ((F0_cands * harmonic >= F_min - bin_tol*fft_bin) & (F0_cands * harmonic <= F_max + bin_tol*fft_bin))

            DM_limit = cand_tools.DM_curve(pm, snr_limit=3)
            dm_offset =  (pm.prop_effect.DM - candidates['dm'])
            dm_cond = np.abs(dm_offset) <= DM_limit

            matched_candidates = candidates[freq_cond & dm_cond]
            matched_candidates = matched_candidates.sort_values(by='snr', key=abs, ascending=False)
            matched_candidates = matched_candidates.head(max_folds)
            matched_candidates['PSR_ID'] = pm.ID
            matched_pulsar_cands.append(matched_candidates)
            print_exe(f'... {len(matched_candidates)} found.')

        return matched_pulsar_cands

    def generate_fold_files(self):
        mode = self.get_mode()

        if len(self.comb_matched_cands) == 0:
            self.write_fold_plan([], f'{self.results_dir}/processing')
            return

        generator = self.fold_file_generators.get(mode)
        if generator is None:
            sys.exit(f"Unknown candidate_matcher_args mode '{mode}'.")
        generator()



if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Search candidate matcher for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    cm_exec = CandidateMatcher(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    cm_exec.setup()
    cm_exec.process()
    cm_exec.generate_fold_files()
