import os
import sys
import glob
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Manager, Pool

import numpy as np
import matplotlib.pyplot as plt

from TOOLS_ar import ARProcessor
from TOOLS_io import parse_cand_file, parse_par_file, parse_JSON, rsync

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from injector.io_tools import merge_filterbanks, FilterbankReader, print_exe


class PulsarxParFolder:
    def __init__(self, tag, processing_args, out_dir, work_dir, mode='INIT'):
        self.processing_args_path = processing_args
        self.processing_args = parse_JSON(processing_args)['nullsar']

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.tag = tag
        self.mode = mode
        self.processing_dir = f'{self.out_dir}/{self.tag}'
        
        manager = Manager()
        self.archive = manager.dict()

    def setup(self):
        self.check_SNR()

        self.create_zap_sting()
        self.transfer_data()

    def check_SNR(self):
        if self.mode != 'INIT':
            init_ar_data =  parse_JSON(f"{self.processing_dir}/02_INIT/INIT_fold_params.json")
            SNR_limit = self.processing_args.get('SNR_limit', 15)
            
            for par_file  in list(self.processing_args['par_files']):
                psr_ID = Path(par_file).stem
                SNR = init_ar_data[psr_ID]['SNR']
                if SNR < SNR_limit:
                    self.processing_args['par_files'].remove(par_file)
                    
        if len(self.processing_args['par_files']) == 0:
            print_exe('No pulsars to null.')
            sys.exit(0)

    def create_zap_sting(self):
        cmask = self.processing_args['fold_pars'].get('channel_mask', '')
        if cmask:
            cmask = cmask.strip()
            self.zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        else:
            self.zap_string = ''

    def transfer_data(self):
        files_dir = f'{self.processing_dir}/01_FILES'

        if self.mode == 'INIT':
            if self.processing_args.get('filtool', False):
                self.data = glob.glob(f'{files_dir}/{self.tag}_FILTOOL*.fil')[0]

            else:
                with open(f'{files_dir}/files.txt') as f:
                    self.data = [line.strip() for line in f if line.strip()]

        elif self.mode == 'OPTIMISE':
            self.data = glob.glob(f'{self.processing_dir}/02_INIT/*INIT.fil')[0]

        elif self.mode == 'CONFIRM':
            self.data = glob.glob(f'{self.processing_dir}/03_OPT/*NULL.fil')[0]

        if type(self.data) == list:
            fb_names = [Path(fb).stem for fb in self.data]

            prefix = os.path.commonprefix(fb_names)
            new_fb_path = f"{self.work_dir}/{prefix}_MERGED_{self.mode}.fil"
            
            self.data.sort()
            merge_filterbanks(self.data, new_fb_path)
            self.data = new_fb_path

    def get_psr_params(self, par_file):
        if self.get_folding_alg(par_file) == '--parfile':
            psr = parse_par_file(par_file)
            p0 = float(psr.get('P0', 0))
            f0 = float(psr.get('F0', 0))
            if not p0 and f0:
                p0 = 1/f0
        else:
            psr = parse_cand_file(par_file)
            p0 = 1/psr['F0']

        blocksize_plan = self.processing_args['fold_pars'].get('blocksize_plan', [2, 10, 0.1])
        input_blocksize = self.processing_args['fold_pars'].get('blocksize', False)
        if not input_blocksize:
            if p0 > blocksize_plan[2]:
                block_size = blocksize_plan[1]
            else:
                block_size = blocksize_plan[0]
        else:
            block_size = input_blocksize

        psrid = Path(par_file).stem
        return psrid, block_size

    def get_folding_alg(self, par_file):
        suffix = Path(par_file).suffix 
        if suffix == '.candfile':
            return '--candfile'
        else:
            return '--parfile'

    def adjust_par_file(self, par_file, psr_id):
        params = parse_par_file(par_file)
        init_ar_data = parse_JSON(f"{self.processing_dir}/02_INIT/INIT_fold_params.json")
        
        params['DM'] = init_ar_data[psr_id]['DM']
        new_par_file = f'{self.work_dir}/{psr_id}_new_parfile.par'
        with open(new_par_file, "w") as f:
            for key, value in params.items():
                f.write(f"{key} {value}\n")

        return new_par_file
    
    def run_parfold(self, par_file):
        fold_args = self.processing_args['fold_pars']

        psr_id, block_size = self.get_psr_params(par_file)
        alg_cmd = self.get_folding_alg(par_file)


        if self.mode == 'INIT':
            # search = '--nof0search --nof1search' 
            search = '--nosearch'
        elif self.mode == 'OPTIMISE':
            search = '--nosearch'
            par_file = self.adjust_par_file(par_file, psr_id)
        else:
            search = ''

        fb = FilterbankReader(self.data, load_fb_stats=(128, 6))
        t_subint = fb.obs_len / fold_args['n_subint']

        tmp_cwd = f'{self.work_dir}/process_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)
        cmd = f"{fold_args['mode']} {search} -o {tmp_cwd}/ --tsubint {t_subint} --nsubband {fb.nchans} -f {self.data} --template {fold_args['template']} {alg_cmd} {par_file} --blocksize {block_size} {self.zap_string} --saveimage"
    
        for flag in fold_args['cmd_flags']:
            if flag != '--saveimage':
                cmd += f" {flag}"

        for key, value in fold_args['cmd'].items():
            if key in ['tsubint', 'nsubband']:
                print(f'{key} parameter not available in Nullsar')
            else:
                cmd += f" --{key} {value}"
        
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

    def run_fold(self, ncpus):
        args = self.processing_args['par_files']

        with Pool(ncpus) as p:
            p.map(self.run_parfold, args)

    def transfer_products(self):

        for par_file in self.processing_args['par_files']:
            psr_id = Path(par_file).stem
            tmp_cwd = f'{self.work_dir}/process_{psr_id}'

            if self.mode == 'INIT':
                files_dir = f'{self.processing_dir}/02_INIT/'
            elif self.mode == 'OPTIMISE':
                files_dir = f'{self.processing_dir}/03_OPT/'
            else:
                files_dir = f'{self.processing_dir}/04_CONFIRM/'

            folds_dir = f'{files_dir}/FOLDS'
            png_path = glob.glob(f'{tmp_cwd}/*.png')
            if png_path:
                rsync(png_path[0], f"{folds_dir}/{psr_id}_mode_{self.mode}.png")

            fits_path = glob.glob(f'{tmp_cwd}/*.px')
            if fits_path:
                rsync(fits_path[0], f"{folds_dir}/{psr_id}_mode_{self.mode}.fits")


        if self.mode == 'CONFIRM':
            files_dir = f'{self.processing_dir}/01_FILES'
            if self.processing_args['delete_filtool']:
                filtool_fb = glob.glob(f'{files_dir}/{self.tag}_FILTOOL*.fil')
                if filtool_fb:
                    os.remove(filtool_fb[0])

            files_dir = f'{self.processing_dir}/02_INIT'
            if self.processing_args['delete_processing']:
                init_fb = glob.glob(f'{files_dir}/*INIT.fil')
                if init_fb:
                    os.remove(init_fb[0])

            for par_file in self.processing_args['par_files']:
                psr_id = Path(par_file).stem
                tmp_cwd = f'{self.work_dir}/process_{psr_id}'
                fits_path = glob.glob(f'{tmp_cwd}/*.px')
                if fits_path:
                    gen_plot(psr_id, self.processing_dir)


def gen_plot(psr_ID, processing_dir):
    fits_path_INIT = f'{processing_dir}/02_INIT/FOLDS/{psr_ID}_mode_INIT.fits'
    fits_path_NULL = f'{processing_dir}/04_CONFIRM/FOLDS/{psr_ID}_mode_CONFIRM.fits'

    archive_INIT = ARProcessor(fits_path_INIT)
    archive_NULL = ARProcessor(fits_path_NULL)

    archives = [archive_INIT, archive_NULL]
    titles = ['INIT', 'NULL']

    fig, axes = plt.subplots(3, 2, figsize=(10, 8), sharex='col')

    for col, (archive, title) in enumerate(zip(archives, titles)):

        IP = archive.get_intensity_prof()
        FP = archive.get_freq_phase()
        TP = archive.get_time_phase()

        if np.max(FP) != 0:
            FP = FP / np.max(FP)

        phase_IP = np.linspace(0, 1, len(IP))

        nchans = FP.shape[0]
        Tobs = TP.shape[0]

        axes[0, col].plot(phase_IP, IP)
        axes[0, col].set_title(f'{title}, S/N: {archive.get_SNR():.2f}')
        axes[0, col].set_ylabel('Intensity')

        axes[1, col].imshow(
            FP,
            origin='lower',
            aspect='auto',
            extent=[0, 1, 0, nchans]
        )
        axes[1, col].set_ylabel('Channel')

        axes[2, col].imshow(
            TP,
            origin='lower',
            aspect='auto',
            extent=[0, 1, 0, Tobs]
        )
        axes[2, col].set_ylabel('Time')
        axes[2, col].set_xlabel('Phase')

    save_path = f'{processing_dir}/04_CONFIRM/CONFIRM_{psr_ID}.png'
    plt.savefig(save_path, dpi=200, bbox_inches='tight')





if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Nullsar folder',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='file tag')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--mode', metavar='int', type=str, required=False, default="INIT", help='folding mode (INIT, OPTIMISE, CONFIRM)')
    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of cpus for injection')

    args = parser.parse_args()
    fold_exec = PulsarxParFolder(args.tag, args.processing_args, args.out_dir, args.work_dir, args.mode)

    fold_exec.setup()
    fold_exec.run_fold(args.ncpus)
    fold_exec.transfer_products()