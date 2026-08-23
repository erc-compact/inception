import glob
import argparse
import subprocess
import numpy as np
from pathlib import Path

import pipeline_tools as inj_tools

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from injector.io_tools import FilterbankReader, print_exe


class PeasoupProcess:
    def __init__(self, process_tag, processing_args, out_dir, work_dir, injection_number):
        self.process_tag = process_tag

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def peasoup_setup(self):
        self.get_injection_report()
        self.parse_tag()
        self.transfer_data()
        self.calc_args()
        self.generate_chan_mask()
        self.generate_birdie_list()
        self.create_dm_list()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        self.report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(self.report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

        self.cdm = self.injection_report['pulsars'][0]['cDM']

    def parse_tag(self):
        tag_split  = self.process_tag.split('_')
        _, _, _, self.tscrunch, _, self.seg_i, self.seg_size = tag_split 
        self.tscrunch = int(self.tscrunch)
        self.seg_i = int(self.seg_i)
        self.seg_size = int(self.seg_size)
        
        tscrunch_list = self.processing_args.get('filtool_args', {'tscrunch': [1]}).get('tscrunch', [1])
        self.tscrunch_index = tscrunch_list.index(self.tscrunch)

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        self.data_inj = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]

        if self.processing_args['peasoup_args']['filtool']:
            data = glob.glob(f"{results_dir}/processing/FILTOOL/*_{self.inj_id}_FILTOOL_0{self.tscrunch_index+1}.fil")[0]
        else:
            data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]
            
        if self.processing_args['peasoup_args'].get('transfer_TMP', True):
            inj_tools.rsync(data, self.work_dir)
            self.data = f'{self.work_dir}/{Path(data).name}'
        else:
            self.data = data
        
    def calc_args(self):
        self.seg_args = self.processing_args['peasoup_args']['segment_plan'][f'{self.seg_size}']

        fb_reader = FilterbankReader(self.data, stats_samples=0)
        fd = self.processing_args.get('filtool_args', {"cmd": {}})['cmd'].get('fd', 1)
        self.default_dedisp_gulp = int((2048.0 / (fb_reader.nchans/fd)) * 1e6)

        seg_fftsize = self.seg_args.get('fftsize', [2, 3, 5])
        if type(seg_fftsize) == list:
            n_sample = fb_reader.n_samples // (self.tscrunch * self.seg_size)
            self.fft_size = inj_tools.next_fast_len(n_sample, seg_fftsize)
        else:
            self.fft_size = seg_fftsize

        ram_limit = self.seg_args.get('ram_limit_gb', 'default')
        if ram_limit == 'default':
            fb_size_GB = self.fft_size * fb_reader.nchans * fb_reader.nbits * 1e-9 
            self.ram_limit_gb = fb_size_GB * 5
        else:
            self.ram_limit_gb = ram_limit

        self.start_sample = int(np.floor(self.seg_i * fb_reader.n_samples / self.seg_size))

    def generate_chan_mask(self):
        chan_mask_string = self.processing_args['peasoup_args'].get('channel_mask', None)
        if chan_mask_string:
            chan_mask_file = f'{self.work_dir}/channel_mask.ascii'

            fb_reader = FilterbankReader(self.data, stats_samples=0)
            ftop = fb_reader.ftop
            fbottom = fb_reader.fbottom
            nchans = fb_reader.nchans

            chan_mask = np.ones(nchans)
            for val in chan_mask_string.split(','):
                if len(val.split(":")) == 1:
                    rstart = float(val)
                    rend = float(val)
                elif len(val.split(":")) == 2:
                    rstart = float(val.split(":")[0])
                    rend = float(val.split(":")[1])
                chbw = (ftop - fbottom) / nchans
                idx0 = int(min(max((rstart - fbottom) // chbw, 0), nchans - 1))
                idx1 = int(max(min(int((rend - fbottom) / chbw + 0.5), nchans - 1), 0))
                chan_mask[idx0:idx1 + 1] = 0

            np.savetxt(chan_mask_file, chan_mask, fmt='%d')
            self.channel_mask =  f' -k {chan_mask_file}'
        else:
            self.channel_mask = ''
    
    def generate_birdie_list(self):
        birdie_list_string = self.processing_args['peasoup_args'].get('birdie_list', None)
        if birdie_list_string:
            birdie_list_file = f'{self.work_dir}/birdie_list.ascii'

            birdies, birdies_width = [], []
            for val in birdie_list_string.split(','):
                f = val.split(":")[0]
                w = val.split(":")[1]
                birdies.append(f)
                birdies_width.append(w)

            np.savetxt(birdie_list_file, np.c_[np.array(birdies, dtype=float),
                                np.array(birdies_width, dtype=float)], fmt="%.2f")
            self.birdie_list = f' -z {birdie_list_file}'
        else:
            self.birdie_list = ''

    def create_dm_list(self):
        self.ddplan = self.processing_args['peasoup_args']['ddplan'][f'{self.tscrunch}']
        
        if self.processing_args['peasoup_args']['inj_DM']:
            DM_values = np.array([pulsar['DM'] for pulsar in self.injection_report['pulsars']])
            if len(DM_values) == 0:
                print_exe('No DMs to search.')
                sys.exit(0)
            else:
                DM_values.sort()
        else:
            low_dm, high_dm, dm_step = self.ddplan
            n_trial = int(round((high_dm-low_dm)/dm_step))
            DM_values = np.linspace(low_dm, high_dm, n_trial, endpoint=False)

        DM_file = f'{self.work_dir}/dm_list_T0{int(self.tscrunch_index)}.ascii'
        np.savetxt(DM_file, DM_values, fmt='%.3f')
        self.DM_file = DM_file
        
    def run_peasoup(self):

        cmd = f"peasoup -i {self.data} --dm_file {self.DM_file} --cdm {self.cdm} --min_snr {self.seg_args['min_snr']} --acc_start {self.seg_args['acc_start']} --acc_end {self.seg_args['acc_end']} -o {self.work_dir} {self.channel_mask} {self.birdie_list}" 

        cmd_args = self.processing_args['peasoup_args']['cmd']
        cmd_args['fft_size'] = self.fft_size
        cmd_args['ram_limit_gb'] = self.ram_limit_gb
        cmd_args['dedisp_gulp'] = cmd_args.get('dedisp_gulp', self.default_dedisp_gulp)

        if not ((self.seg_i == 0) and (self.seg_size == 1)):
            cmd_args['start_sample'] = self.start_sample

        for key, value in cmd_args.items():
            if key in ['end_sample']:
                pass
            cmd += f" --{key} {value}"
              
        subprocess.run(cmd, shell=True)
    
    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        processing_dir = f'{results_dir}/processing'

        prefix = f"{self.processing_args['injection_args']['id']}_{self.process_tag}"

        xml_name_old = f'{self.work_dir}/overview.xml'
        if os.path.exists(xml_name_old):
            if self.processing_args['peasoup_args']['save_xml']:
                inj_tools.rsync(xml_name_old, f"{processing_dir}/PEASOUP/{prefix}.xml")
        else:
            print_exe('No xml file produced!')

            
        peasoup_cands = f'{self.work_dir}/candidates.peasoup'
        if os.path.exists(peasoup_cands):
            if self.processing_args['peasoup_args']['save_peasoup']:
                inj_tools.rsync(peasoup_cands, f'{processing_dir}/PEASOUP/{prefix}.peasoup')
        else:
            print_exe('No candidates.peasoup produced!')

        if self.processing_args['peasoup_args']['save_birdie']:
            inj_tools.rsync(f'{self.work_dir}/birdie_list.ascii', f'{processing_dir}/PEASOUP/{prefix}_birdie.ascii')

        if self.processing_args['peasoup_args']['save_cmask']:
            inj_tools.rsync(f'{self.work_dir}/channel_mask.ascii', f'{processing_dir}/PEASOUP/{prefix}_cmask.ascii')

        if self.processing_args['peasoup_args']['save_dm_list']:
            inj_tools.rsync(self.DM_file, f'{processing_dir}/PEASOUP/{prefix}_{Path(self.DM_file).name}')

        if self.processing_args['peasoup_args']['delete_filtool_fb']:
            data = glob.glob(f"{processing_dir}/*_{self.inj_id}_FILTOOL_0{self.tscrunch_index+1}.fil")
            if data:
                os.remove(data[0])

        if self.processing_args['peasoup_args']['delete_inj_fb']:
            data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")
            if data:
                os.remove(data[0])

        


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='peasoup for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--tag', metavar='str', required=True, type=str, help='search process tag')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    peasoup_exec = PeasoupProcess(args.tag, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    peasoup_exec.peasoup_setup()
    peasoup_exec.run_peasoup()
    peasoup_exec.transfer_products()