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
    def __init__(self, tscrunch_index, processing_args, out_dir, work_dir, injection_number):
        self.tscrunch_index = tscrunch_index

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def peasoup_setup(self):
        self.get_injection_report()
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

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        self.data_inj = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]

        if self.processing_args['peasoup_args']['filtool']:
            data = glob.glob(f"{results_dir}/processing/*_{self.inj_id}_FILTOOL_0{self.tscrunch_index+1}.fil")[0]
            inj_tools.rsync(data, self.work_dir)

        else:
            data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]
            inj_tools.rsync(data, self.work_dir)

        self.data = f'{self.work_dir}/{Path(data).name}'
        

    def calc_args(self):
        fb_reader = FilterbankReader(self.data_inj, stats_samples=0)
        fd = self.processing_args.get('filtool_args', {"cmd": {}})['cmd'].get('fd', 1)
        self.default_dedisp_gulp = int((2048.0 / (fb_reader.nchans/fd)) * 1e6)

        self.default_fft_size = 2**int(np.ceil(np.log2(fb_reader.n_samples)))

        fb_size_GB = self.default_fft_size * fb_reader.nchans * fb_reader.nbits * 1e-9 
        self.default_ram_limit_gb = fb_size_GB * 3

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
        ddplan_list = inj_tools.create_DDplan(self.processing_args['peasoup_args']['ddplan'])
        self.ddplan = ddplan_list[self.tscrunch_index]
        self.tscrunch = int(self.ddplan.tscrunch)
        
        if self.processing_args['peasoup_args']['inj_DM']:
            DM_values = [pulsar['DM'] for pulsar in self.injection_report['pulsars']]
            DM_values = DM_values[(DM_values >= self.ddplan.low_dm) & (DM_values <= self.ddplan.high_dm)]
            if len(DM_values) == 0:
                print_exe('No DMs to search.')
                sys.exit(0)
            else:
                DM_values.sort()
        else:
            n_trial = int(round((self.ddplan.high_dm-self.ddplan.low_dm)/self.ddplan.dm_step))
            if self.tscrunch_index == len(ddplan_list)-1:
                n_trial, endpoint = (n_trial + 1, True)
            else:
                n_trial, endpoint = (n_trial, False)

            DM_values = np.linspace(self.ddplan.low_dm, self.ddplan.high_dm, n_trial, endpoint=endpoint)

        DM_file = f'{self.work_dir}/dm_list_T0{int(self.tscrunch_index)}.ascii'
        np.savetxt(DM_file, DM_values, fmt='%.3f')
        self.DM_file = DM_file
        
    def run_peasoup(self):

        cmd = f"peasoup -i {self.data} --dm_file {self.DM_file} -o {self.work_dir} {self.channel_mask} {self.birdie_list}" 

        cmd_args = self.processing_args['peasoup_args']['cmd']
        cmd_args['fft_size'] = int(cmd_args.get('fft_size', self.default_fft_size) // self.tscrunch)
        cmd_args['dedisp_gulp'] = cmd_args.get('dedisp_gulp', self.default_dedisp_gulp)
        cmd_args['ram_limit_gb'] = cmd_args.get('ram_limit_gb', self.default_ram_limit_gb)

        for key, value in cmd_args.items():
            cmd += f" --{key} {value}"
              
        subprocess.run(cmd, shell=True)

    def process_candidates(self, processing_dir, xml_name_old, prefix):
        save_csv = self.processing_args['peasoup_args']['save_csv']
        match_inj =self.processing_args['peasoup_args']['candidate_matcher']
        if match_inj['match_inj'] or save_csv:
            from candidate_tools import xml2csv

            csv_cands = xml2csv(xml_name_old, f"{processing_dir}/{prefix}.csv", save_csv)

        if match_inj['match_inj']:
            from candidate_tools import CandMatcher

            fft_size = int(self.processing_args['peasoup_args']['cmd'].get('fft_size', self.default_fft_size // self.tscrunch) * self.tscrunch)
            ephem = self.processing_args['injection_args']['ephem']
            if ephem != 'builtin':
                inj_tools.rsync(ephem, self.work_dir)
                ephem = f'./{Path(ephem).name}'
            cand_matcher = CandMatcher(self.report_path, csv_cands, self.data_inj, fft_size, ephem, corr_period=True)

            candidate_root = f"{processing_dir}/{self.processing_args['injection_args']['id']}_{self.inj_id}_{match_inj['tag']}_0{self.tscrunch_index+1}"
            cand_matcher.generate_files(candidate_root, max_cand_per_inj=match_inj.get('n_cands_per_inj', 1), 
                                        pepoch_ref=0.5, snr_limit=match_inj.get('DM_snr_limit', 3), max_harmonic=match_inj.get('max_harmonic', 4),
                                        create_candfile=match_inj.get('create_candfile', True))
    

    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        processing_dir = f'{results_dir}/processing'

        prefix = f"{self.processing_args['injection_args']['id']}_{self.inj_id}_DM_{self.ddplan.low_dm:.6f}_{self.ddplan.high_dm:.6f}"

        xml_name_old = f'{self.work_dir}/overview.xml'
        if os.path.exists(xml_name_old):
            if self.processing_args['peasoup_args']['save_xml']:
                inj_tools.rsync(xml_name_old, f"{processing_dir}/{prefix}.xml")

            self.process_candidates(processing_dir, xml_name_old, prefix)
        else:
            print_exe('No xml file produced!')

            
        peasoup_cands = f'{self.work_dir}/candidates.peasoup'
        if os.path.exists(peasoup_cands):
            if self.processing_args['peasoup_args']['save_peasoup']:
                inj_tools.rsync(peasoup_cands, f'{processing_dir}/{prefix}.peasoup')
        else:
            print_exe('No candidates.peasoup produced!')

        if self.processing_args['peasoup_args']['save_birdie']:
            inj_tools.rsync(f'{self.work_dir}/birdie_list.ascii', f'{processing_dir}/{prefix}_birdie.ascii')

        if self.processing_args['peasoup_args']['save_cmask']:
            inj_tools.rsync(f'{self.work_dir}/channel_mask.ascii', f'{processing_dir}/{prefix}_cmask.ascii')

        if self.processing_args['peasoup_args']['save_dm_list']:
            inj_tools.rsync(self.DM_file, f'{processing_dir}/{prefix}_{Path(self.DM_file).name}')

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
    parser.add_argument('--tscrunch_index', metavar='int', required=False, default=0, type=int, help='tscrunch index of filterbank file to search')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()

    peasoup_exec = PeasoupProcess(args.tscrunch_index, args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    peasoup_exec.peasoup_setup()
    peasoup_exec.run_peasoup()
    peasoup_exec.transfer_products()