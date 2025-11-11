import os
import sys
import glob
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Manager, Pool

from nullsar_tools import parse_cand_file, parse_par_file, parse_JSON, rsync

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
        self.processing_dir = f'{self.out_dir}/PROCESSING/{self.tag}'
        
        manager = Manager()
        self.archive = manager.dict()

    def setup(self):
        self.check_SNR()

        self.create_zap_sting()
        self.transfer_data()

    def check_SNR(self):
        if self.mode != 'INIT':
            files_dir = f'{self.processing_dir}/01_FILES/NULLSAR'
            init_ar_data =  parse_JSON(f"{files_dir}/INIT_fold_params.json")
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
                self.data = glob.glob(f'{files_dir}/*FILTOOL*.fil')[0]

            else:
                with open(f'{files_dir}/files.txt') as f:
                    self.data = [line.strip() for line in f if line.strip()]

        elif self.mode == 'OPTIMISE':
            self.data = glob.glob(f'{files_dir}/NULLSAR/*INIT.fil')[0]

        elif self.mode == 'CONFIRM':
            self.data = glob.glob(f'{files_dir}/NULLSAR/*NULL.fil')[0]

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
    
    def run_parfold(self, par_file):
        fold_args = self.processing_args['fold_pars']

        psr_id, block_size = self.get_psr_params(par_file)
        alg_cmd = self.get_folding_alg(par_file)

        search = '--nosearch' if (self.mode != 'CONFIRM') else ''

        fb = FilterbankReader(self.data)
        t_subint = fb.obs_len / fold_args['n_subint']

        tmp_cwd = f'{self.work_dir}/process_{psr_id}'
        os.makedirs(tmp_cwd, exist_ok=True)
        cmd = f"{fold_args['mode']} {search} -o {tmp_cwd}/ --t_subint {t_subint} --nsubband {fb.nchans} -f {self.data} --template {fold_args['template']} {alg_cmd} {par_file} --blocksize {block_size} {self.zap_string} --saveimage"
    
        for flag in fold_args['cmd_flags']:
            if flag != '--saveimage':
                cmd += f" {flag}"

        for key, value in fold_args['cmd'].items():
            if key in ['t_subint', 'nsubband']:
                print(f'{key} parameter not available in Nullsar')
            else:
                cmd += f" --{key} {value}"
        
        subprocess.run(cmd, shell=True, cwd=tmp_cwd)

    def run_fold(self, ncpus):
        args = self.processing_args['par_files']

        with Pool(ncpus) as p:
            p.map(self.run_parfold, args)

    def transfer_products(self):
        files_dir = f'{self.processing_dir}/01_FILES/NULLSAR/'
        folds_dir = f'{files_dir}/FOLDS'
        os.makedirs(folds_dir, exist_ok=True)

        for par_file in self.processing_args['par_files']:
            psr_id = Path(par_file).stem
            tmp_cwd = f'{self.work_dir}/process_{psr_id}'
            png_path = glob.glob(f'{tmp_cwd}/*.png')
            if png_path:
                rsync(png_path[0], f"{folds_dir}/{psr_id}_mode_{self.mode}.fits")

            if self.mode != 'CONFIRM':
                fits_path = glob.glob(f'{tmp_cwd}/*.px')
                rsync(fits_path[0], f"{folds_dir}/{psr_id}_mode_{self.mode}.fits")

        if self.mode == 'CONFIRM':
            if self.processing_args['delete_filtool']:
                filtool_fb = glob.glob(f'{files_dir}/*FILTOOL*.fil')
                if filtool_fb:
                    os.remove(filtool_fb[0])

            if self.processing_args['delete_processing']:
                init_fb = glob.glob(f'{files_dir}/*INIT.fil')
                if init_fb:
                    os.remove(init_fb[0])

                fits_files = glob.glob(f'{folds_dir}/*.fits')
                for file in fits_files:
                    os.remove(file)
            


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