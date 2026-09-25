import os
import sys
import glob
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'PY_general')))
import pipeline_tools as inj_tools

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from injector.io_tools import FilterbankReader, print_exe


class PrestoFoldCandProcess:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.injection_number = injection_number
        self.results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

    def fold_setup(self):
        self.get_injection_report()
        self.get_candidates()
        self.mode = self.processing_args['presto_candfold_args'].get('fold_mode', 'filterbank')
        if self.mode != 'dat':
            self.transfer_filterbank()

    def get_injection_report(self):
        report_path = glob.glob(f'{self.results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def get_candidates(self):
        processing_dir = f'{self.results_dir}/processing/PRESTO'
        cand_csv = glob.glob(f'{processing_dir}/{self.process_tag}.csv')
        if cand_csv:
            self.candidates = pd.read_csv(cand_csv[0], index_col=0)
        else:
            print_exe('No matched candidate file found.')
            sys.exit(0)

    def transfer_filterbank(self):
        data = glob.glob(f"{self.results_dir}/*_{self.inj_id}.fil")
        if not data:
            print_exe('No injected filterbank found.')
            sys.exit(0)

        if self.processing_args['presto_candfold_args'].get('transfer_TMP', True):
            inj_tools.rsync(data[0], self.work_dir)
            self.data = f'{self.work_dir}/{Path(data[0]).name}'
        else:
            self.data = data[0]

    def get_mask(self):
        f_args = self.processing_args['presto_candfold_args']
        presto_out_dir = f'{self.results_dir}/processing/PRESTO'
        if f_args['mask'] == 'rfifind':
            return f'-mask {presto_out_dir}/{self.inj_id}_rfifind.mask'
        elif f_args['mask']:
            return f"-mask {f_args['mask']}"
        else:
            return ''

    @staticmethod
    def fdd_flag(cand):
        return f"-fdd {cand['F2']}" if cand['F2'] else ''

    def get_segment_bounds(self, seg_i, seg_n):
        fb_reader = FilterbankReader(self.data, stats_samples=0)
        n_samples = fb_reader.n_samples

        start_sample = int(np.floor(seg_i * n_samples / seg_n))
        end_sample = int(np.floor((seg_i + 1) * n_samples / seg_n))

        start = max(start_sample / n_samples, 0)
        end = min(end_sample / n_samples, 1)
        return start, end

    def bary_flag(self):
        """prepfold has no -nobary (that is prepdata's spelling); it opts out of
        barycentring with -topo. Defaults to whatever the search used, since the
        candidate's F0/F1/F2 are in the frame the time series was dedispersed in."""
        f_args = self.processing_args['presto_candfold_args']
        search_bary = self.processing_args['presto_search_args'].get('bary', False)
        return '' if f_args.get('bary', search_bary) else '-topo'

    def fold_candidate_fb(self, cand):
        f_args = self.processing_args['presto_candfold_args']
        bary = self.bary_flag()
        mask = self.get_mask()

        start, end = self.get_segment_bounds(int(cand['seg_i']), int(cand['seg_n']))

        cwd = f"{self.work_dir}/_{cand['PSR_ID']}"
        os.makedirs(cwd, exist_ok=True)
        out_file = f"{cwd}/{cand['PSR_ID']}_{cand.name}"

        cmd = (f"prepfold {bary} -noxwin -o {out_file} -f {cand['F0']} -fd {cand['F1']} {self.fdd_flag(cand)} "
               f"-dm {cand['dm']} -start {start} -end {end} {mask} {self.data}")

        cmd = inj_tools.add_cmd_args(cmd, f_args, skip_flags=['-noxwin', '-topo'],
                                     skip_keys=['o', 'f', 'fd', 'fdd', 'dm', 'start', 'end', 'mask'])

        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=cwd)

    def fold_candidate_dat(self, cand):
        dat_dir = f'{self.results_dir}/processing/PRESTO/FOLD_DAT'
        root = f"{self.inj_id}_SEG_{int(cand['seg_i'])}_{int(cand['seg_n'])}_DS{int(cand['downsample'])}_DM{cand['dm']:.2f}"

        dat_file = glob.glob(f'{dat_dir}/{root}.dat')
        inf_file = glob.glob(f'{dat_dir}/{root}.inf')
        if not (dat_file and inf_file):
            print_exe(f'No saved segment .dat found for {root}, skipping {cand["PSR_ID"]}.')
            return

        cwd = f"{self.work_dir}/_{cand['PSR_ID']}"
        os.makedirs(cwd, exist_ok=True)
        inj_tools.rsync(dat_file[0], cwd)
        inj_tools.rsync(inf_file[0], cwd)

        out_file = f"{cwd}/{cand['PSR_ID']}_{cand.name}"
        # no bary flag here: the .dat was already barycentred (or not) by prepdata
        # and its .inf records which, so prepfold must not be told to re-decide
        cmd = (f"prepfold -noxwin -o {out_file} -f {cand['F0']} -fd {cand['F1']} {self.fdd_flag(cand)} "
               f"-dm {cand['dm']} {cwd}/{root}.dat")

        cmd = inj_tools.add_cmd_args(cmd, self.processing_args['presto_candfold_args'],
                                     skip_flags=['-noxwin', '-topo'],
                                     skip_keys=['o', 'f', 'fd', 'fdd', 'dm', 'start', 'end', 'mask'])

        inj_tools.print_exe(cmd)
        subprocess.run(cmd, shell=True, cwd=cwd)

    def fold_candidate(self, cand):
        if self.mode == 'dat':
            self.fold_candidate_dat(cand)
        else:
            self.fold_candidate_fb(cand)

    def run_fold(self, ncpus):
        cands = [row for _, row in self.candidates.iterrows()]
        with Pool(ncpus) as p:
            p.map(self.fold_candidate, cands)

    # prepfold appends its own '_<period>ms_Cand' suffix to -o, so every product is
    # renamed to a deterministic '{PSR_ID}_CAND{index}_...' name the collector can
    # glob for. Longest suffix first: '.pfd.bestprof' must win over '.pfd'.
    CAND_PRODUCTS = [('.pfd.bestprof', 'save_bestprof', '.bestprof'),
                     ('.pfd.ps',       'save_ps',       '.ps'),
                     ('.pfd',          'save_pfd',      '.pfd'),
                     ('.png',          'save_png',      '.png')]

    def transfer_products(self):
        results_dir = f'{self.results_dir}/inj_cands/PRESTO/{self.process_tag}'
        os.makedirs(results_dir, exist_ok=True)

        f_args = self.processing_args['presto_candfold_args']
        tag = f"{self.processing_args['injection_args']['id']}_{self.inj_id}_inj_{self.injection_number:06}"

        for idx, cand in self.candidates.iterrows():
            psr_id = cand['PSR_ID']
            cwd = f'{self.work_dir}/_{psr_id}'
            if not os.path.isdir(cwd):
                continue

            prefix = f'{psr_id}_{idx}_'
            for filename in os.listdir(cwd):
                if not filename.startswith(prefix):
                    continue

                for suffix, save_flag, out_ext in self.CAND_PRODUCTS:
                    if not filename.endswith(suffix):
                        continue
                    if f_args.get(save_flag, True):
                        shutil.move(f'{cwd}/{filename}',
                                    f'{results_dir}/{psr_id}_CAND{idx}_{tag}{out_ext}')
                    break


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Presto candidate-folder for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--ncpus', metavar='int', type=int, required=False, default=1, help='number of cpus to use')
    args = parser.parse_args()

    fold_exec = PrestoFoldCandProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    fold_exec.fold_setup()
    fold_exec.run_fold(args.ncpus)
    fold_exec.transfer_products()
