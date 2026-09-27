import re
import glob
import json
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
    # PSR_ID for candidates folded under fold_all that matched no injected pulsar
    UNMATCHED = 'UNMATCHED'

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

    def get_data_frame(self):
        """Frame of the time series the search actually ran on - not the frame the
        pulsar's spin model is defined in. peasoup searches the raw filterbank, so
        Earth motion is still present; PRESTO's prepdata removes it unless
        presto_search_args['bary'] is false."""
        if self.get_mode() == 'PRESTO':
            return 'bary' if self.processing_args['presto_search_args'].get('bary', False) else 'topo'

        return 'topo'


    @staticmethod
    def parse_peasoup_tag(xml_tag):
        splits = xml_tag.split('_')
        return splits[-2], splits[-1], splits[-4]

    def correct_freq(self, csv_cands, fftsize, dt):
        sys.exit(1)
        # cand_tools.correct_fftsize_offset(period, acc, fftsize, nsamples, dt)

    def load_peasoup_candidates(self):
        processing_dir = f'{self.results_dir}/processing'

        xmls = glob.glob(f'{processing_dir}/PEASOUP/*.xml')
        cands_list = []
        for xml in xmls:
            csv_cands, pepoch, fftsize, dt = cand_tools.xml2csv(xml)
            if pepoch is not None:
                csv_cands['pepoch'] = pepoch
            else:
                self.correct_freq(csv_cands, fftsize, dt)
            csv_cands['n_samples'] = fftsize
            csv_cands['T'] = fftsize * dt
            csv_cands['XML'] = Path(xml).name
            s0, s1, tscrunch = self.parse_peasoup_tag(Path(xml).stem)
            csv_cands['tscrunch'] = tscrunch
            csv_cands['downsample'] = int(tscrunch)
            csv_cands['segment'] = f'{s0}_{s1}'
            # how PulsarX's folds of these are named: the candidate id in this XML,
            # plus its DDPLAN, since one MATCHED_SEG folder holds every downsample
            csv_cands['fold_id'] = f'DDPLAN_{tscrunch}_' + csv_cands['xml_id'].astype(str)
            csv_cands['F_match'] = 1 / csv_cands['period']
            cands_list.append(csv_cands)

        candidates = pd.concat(cands_list, ignore_index=True) if cands_list else pd.DataFrame()

        if len(candidates):
            candidates.to_csv(f'{processing_dir}/PEASOUP/inj_{self.injection_number:06}_PEASOUP_candidates.csv')
        return candidates

    def generate_peasoup_fold_files(self):
        processing_dir = f'{self.results_dir}/processing/PEASOUP'
        unique_segments = np.unique(self.fold_cands['segment'])

        process_tags = []
        for segment in unique_segments:
            s0, s1 = segment.split('_')
            segment_cands = self.fold_cands[self.fold_cands['segment'] == segment]
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
        if len(candidates) == 0:
            print_exe('PRESTO found no candidates.')
            return pd.DataFrame()

        # presto_candfold names its folds by the row in the candidates csv
        candidates['fold_id'] = 'CAND' + candidates.index.astype(str)

        ref_obs = self.setup_manager.pulsar_models[0].obs
        mid_time_sec = (candidates['seg_i'] + 0.5) / candidates['seg_n'] * ref_obs.obs_len
        candidates['pepoch'] = ref_obs.sec2mjd(mid_time_sec.values)

        candidates.to_csv(f'{processing_dir}/inj_{self.injection_number:06}_PRESTO_candidates.csv')
        return candidates

    def generate_presto_fold_files(self):
        processing_dir = f'{self.results_dir}/processing/PRESTO'
        unique_segments = np.unique(self.fold_cands['segment'])

        process_tags = []
        for segment in unique_segments:
            segment_cands = self.fold_cands[self.fold_cands['segment'] == segment]

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

        self.candidates = loader()
        self.dm_windows = {}

        if len(self.candidates) == 0:
            self.matches = pd.DataFrame(columns=['PSR_ID', 'segment'])
            self.fold_cands = pd.DataFrame()
            return

        self.matches = self.match_candidates(self.candidates, harmonics)
        self.fold_cands = self.select_folds(self.candidates, self.matches)

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
        matcher_args = self.processing_args['candidate_matcher_args']
        dm_level = matcher_args.get('dm_level', 0.3)
        data_frame = self.get_data_frame()
        h_arr = np.asarray(harmonics, dtype=float)
        n_cands = len(candidates)
        obs_rv, topo_sec, bary_sec = self.get_obs_params()
        bin_tol = 2

        seg = candidates['segment'].str.split('_', expand=True).astype(int)
        seg_groups = candidates.groupby([seg[0], seg[1]]).indices
        obs_frac = topo_sec / topo_sec[-1]

        print_exe(f'Matching {self.get_mode()} candidates from a {data_frame}centric time series ...')

        fft_bin = 1 / candidates['T'].values
        F0_cands = candidates['F_match'].values

        for pm in self.setup_manager.pulsar_models:
            print_exe(f'Matching PSR {pm.ID} ...')

            frame = pm.pulsar_pars['frame']
            time = bary_sec if frame == 'bary' else topo_sec
            if frame == data_frame:
                psr_rv = np.zeros_like(obs_rv)
            elif frame == 'bary':
                psr_rv = obs_rv.copy()
            else:
                psr_rv = -obs_rv

            rv_PSR = cand_tools.add_PSR_rv_curve(pm, time, psr_rv, pm.pulsar_pars['ACCEPOCH']) 

            F_min = np.empty(len(candidates))
            F_max = np.empty(len(candidates))
            for (seg_i, seg_n), idx in seg_groups.items():
                window = (obs_frac >= seg_i / seg_n) & (obs_frac <= (seg_i + 1) / seg_n)
                F_min[idx], F_max[idx] = cand_tools.get_freq_bounds(rv_PSR[window], pm)

            scaled = np.outer(F0_cands, h_arr)
            in_band = ((scaled >= (F_min - bin_tol * fft_bin)[:, None]) &
                       (scaled <= (F_max + bin_tol * fft_bin)[:, None]))
            freq_cond = in_band.any(axis=1)

            centre = 0.5 * (F_min + F_max)
            best_h = np.argmin(np.where(in_band, np.abs(scaled - centre[:, None]), np.inf), axis=1)
            dF0_bins = (scaled[np.arange(n_cands), best_h] - centre) * candidates['T'].values

            width_low, width_high = cand_tools.dm_match_bounds(pm, dm_level)
            self.dm_windows[pm.ID] = (-width_low, width_high)
            dm_offset = candidates['dm'].values - pm.prop_effect.DM
            dm_cond = (dm_offset >= -width_low) & (dm_offset <= width_high)

            keep = freq_cond & dm_cond
            matched_candidates = candidates[keep].copy()
            matched_candidates['PSR_ID'] = pm.ID
            matched_candidates['harmonic'] = 1.0 / h_arr[best_h[keep]]
            matched_candidates['dF0_bins'] = dF0_bins[keep]
            matched_candidates['dDM'] = dm_offset[keep]
            matched_pulsar_cands.append(matched_candidates)
            print_exe(f'... {len(matched_candidates)} found.')

        return pd.concat(matched_pulsar_cands)

    def select_folds(self, candidates, matches):
        matcher_args = self.processing_args['candidate_matcher_args']
        max_folds = matcher_args['max_folds']
        by_snr = dict(by='snr', key=abs, ascending=False)

        if not matcher_args.get('fold_all', False):
            return (matches.sort_values(**by_snr)
                           .groupby('PSR_ID', sort=False).head(max_folds))

        top = candidates.sort_values(**by_snr).groupby('segment', sort=False).head(max_folds)

        match_cols = ['PSR_ID', 'harmonic', 'dF0_bins', 'dDM']
        matched = top.join(matches[match_cols], how='inner')
        unmatched = top[~top.index.isin(matches.index)].assign(PSR_ID=self.UNMATCHED)

        return pd.concat([matched, unmatched]).sort_values(**by_snr)

    def generate_fold_files(self):
        mode = self.get_mode()

        if len(self.fold_cands) == 0:
            self.write_fold_plan([], f'{self.results_dir}/processing')
            return

        generator = self.fold_file_generators.get(mode)
        if generator is None:
            sys.exit(f"Unknown candidate_matcher_args mode '{mode}'.")
        generator()

    def write_matches(self):
        mode = self.get_mode()
        processing_dir = f'{self.results_dir}/processing/{mode}'
        os.makedirs(processing_dir, exist_ok=True)
        self.matches.to_csv(f'{processing_dir}/inj_{self.injection_number:06}_{mode}_matches.csv')

    def get_segments(self):
        plan_section = {'PEASOUP': 'peasoup_args', 'PRESTO': 'presto_search_args'}.get(self.get_mode())
        plan = self.processing_args.get(plan_section, {}).get('segment_plan', {})

        segments = [f'{si}_{s_plan}' for s_plan in plan for si in range(int(s_plan))]
        if len(self.candidates):
            segments += [s for s in self.candidates['segment'].unique() if s not in segments]
        return segments

    def write_match_summary(self):
        mode = self.get_mode()
        matcher_args = self.processing_args['candidate_matcher_args']
        inj_tag = f'inj_{self.injection_number:06}'
        seg_key = lambda s: f'SEG_{s}'
        natural = lambda s: [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]
        fold_names = {'PRESTO': '<PSR_ID>_<fold_id>_... (<fold_id> = CAND<cand>)',
                      'PEASOUP': '..._<fold_id>.png / .ar (<fold_id> = DDPLAN_<ds>_<id>, <id> being the '
                                 'candidate id in that segment\'s XML ..._DDPLAN_<ds>_SEG_<i>_<n>.xml)'}

        pulsars = []
        for pm in self.setup_manager.pulsar_models:
            psr_matches = self.matches[self.matches['PSR_ID'] == pm.ID]
            psr_folds = self.fold_cands[self.fold_cands['PSR_ID'] == pm.ID] if len(self.fold_cands) else psr_matches[:0]

            folds = {seg_key(s): sorted(set(psr_folds.loc[psr_folds['segment'] == s, 'fold_id']), key=natural)
                     for s in self.get_segments() if (psr_folds['segment'] == s).any()}

            best = None
            if len(psr_matches):
                b_id, b = next(psr_matches.sort_values(by='snr', key=abs, ascending=False).iterrows())
                best = {
                    'cand': int(b_id),
                    'fold_id': b['fold_id'],
                    'segment': seg_key(b['segment']),
                    'snr': round(float(b['snr']), 2),
                    'P0': round(float(b['period']), 9),
                    'DM': round(float(b['dm']), 3),
                    'downsample': int(b['downsample']) if 'downsample' in b else None,
                    'harmonic': float(b['harmonic']),
                    'dF0_bins': round(float(b['dF0_bins']), 2),
                    'dDM': round(float(b['dDM']), 3),
                    'folded': int(b_id) in psr_folds.index,
                }

            dm_window = self.dm_windows.get(pm.ID)
            pulsars.append({
                'PSR_ID': pm.ID,
                'injected': {'P0': round(float(pm.PX_list[0]), 9),
                             'DM': round(float(pm.prop_effect.DM), 3),
                             'SNR': float(pm.SNR)},
                'detected': best is not None,
                'dm_window': [round(float(w), 3) for w in dm_window] if dm_window else None,
                'best': best,
                'n_hits': {seg_key(s): int((psr_matches['segment'] == s).sum()) for s in self.get_segments()},
                'folds': folds,
            })

        summary = {
            'injection': inj_tag,
            'search': mode,
            'fold_all': bool(matcher_args.get('fold_all', False)),
            'max_folds': matcher_args['max_folds'],
            'legend': {
                'cand': f'row in processing/{mode}/{inj_tag}_{mode}_candidates.csv',
                'fold_id': f'how the fold appears in the fold file names: {fold_names.get(mode, "")}; '
                           f'folds lists these',
                'segment': f'SEG_<i>_<n> = segment i (from 0) of n equal parts of the observation; '
                           f'its folds are in inj_cands/{mode}/MATCHED_SEG_<i>_<n>/',
                'harmonic': 'detected frequency / pulsar frequency (2 = found at the 2nd harmonic)',
                'dF0_bins': 'detected minus expected frequency in that segment, in Fourier bins of the search',
                'dDM': 'detected minus injected DM; dm_window is the accepted range',
            },
            'pulsars': pulsars,
        }

        out_dir = f'{self.results_dir}/inj_cands/{mode}'
        os.makedirs(out_dir, exist_ok=True)
        with open(f'{out_dir}/{inj_tag}_{mode}_matches.json', 'w') as f:
            json.dump(summary, f, indent=2)



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
    cm_exec.write_matches()
    cm_exec.write_match_summary()
