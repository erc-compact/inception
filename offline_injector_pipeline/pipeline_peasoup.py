import os
import re
import sys
import argparse
import subprocess
import numpy as np
from pathlib import Path

from pipeline_tools import PipelineTools
from pipeline_inj_cand_sifter import CandFinder
from inception.injector.io_tools import FilterbankReader


class PeasoupExec(PipelineTools):
    def __init__(self, tscrunch_index, search_args, injection_file, out_dir, data_dir, injection_number, n_nearest=-1):
        super().__init__(search_args)
        self.work_dir = os.getcwd()
        self.tscrunch_index = int(tscrunch_index)
        self.tscrunch = self.get_tscrunch()
        self.out_dir = out_dir
        self.data_dir = data_dir
        self.injection_number = injection_number
        self.n_nearest = n_nearest

        self.inj_file = injection_file
        self.fb = self.get_inputs()
        
        self.gulp_size = self.get_gulp_size()

    def get_tscrunch(self):
        DD_plan = self.create_DDplan()
        return [dm_range.tscrunch for dm_range in DD_plan][self.tscrunch_index]
        
    def get_inputs(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        for filename in os.listdir(results_dir):
            if 'report' in filename:
                injection_report = f'{results_dir}/{filename}'

        process_dir = f'{results_dir}/processing'
        filterbank = ''
        for filename in os.listdir(process_dir):
            if filename.endswith(f'{self.tscrunch_index+1}.fil'):
                subprocess.run(f"rsync -Pav {process_dir}/{filename} {self.work_dir}", shell=True)
                filterbank = filename

        self.inj_report_path = injection_report
        self.inj_report = self.parse_JSON(injection_report)
        
        if not filterbank:
            DD_plan = self.create_DDplan()
            xml_name = [f'overview_dm_{dm_range.low_dm:.6f}_{dm_range.high_dm:.6f}.xml' for dm_range in DD_plan][self.tscrunch_index]
            inj_ID = self.inj_report['injection_report']['ID']
            xml_name_new = f'{self.data_ID}_{inj_ID}_{xml_name}'

            self.peasoup_failed(xml_name_new, xml_name, process_dir)
            sys.exit(0)
        
        return filterbank
    
    def get_gulp_size(self):
        fscrunch = int(self.processing_args.get('fscrunch', 1))

        fb_reader = FilterbankReader(self.fb)
        default_gulp_size = int((2048.0 / (fb_reader.nchans / fscrunch)) * 1e6)
        del fb_reader
        return self.processing_args.get('gulp_size', default_gulp_size)
    
    def generate_chan_mask(self, chan_mask_csv, outfile):
        fb_reader = FilterbankReader(self.fb)
        ftop = fb_reader.ftop
        fbottom = fb_reader.fbottom
        nchans = fb_reader.nchans
        del fb_reader
        chan_mask = np.ones(nchans)
        for val in chan_mask_csv.split(','):
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
        np.savetxt(outfile, chan_mask, fmt='%d')
    
    @staticmethod
    def generate_birdie_list(birdie_csv, outfile):
        birdies = []
        birdies_width = []

        for val in birdie_csv.split(','):
            f = val.split(":")[0]
            w = val.split(":")[1]
            birdies.append(f)
            birdies_width.append(w)

        np.savetxt(outfile, np.c_[np.array(birdies, dtype=float),
                               np.array(birdies_width, dtype=float)], fmt="%.2f")
       
    def generate_files(self):
        chan_mask_csv = self.processing_args['channel_mask']
        chan_mask_file = os.path.join(self.work_dir, 'channel_mask.ascii')
        self.generate_chan_mask(chan_mask_csv, chan_mask_file)

        birdie_list_csv = self.processing_args['birdie_list']
        birdie_list_file = os.path.join(self.work_dir, 'birdie_list.ascii')
        self.generate_birdie_list(birdie_list_csv, birdie_list_file)
        return chan_mask_file, birdie_list_file
    
    def find_n_nearest(self, dm_arr, tscrunch_arr, target, n):
        if n == -1:
            dm_nearest, tscrunch_nearest =  dm_arr, tscrunch_arr
        else:
            nearest_indices = np.argsort(np.abs(dm_arr - target))[:n]
            dm_values = dm_arr[nearest_indices]
            tscrunch_values = tscrunch_arr[nearest_indices]
            dm_nearest, tscrunch_nearest = dm_values, tscrunch_values
        
        return dm_nearest[np.where(tscrunch_nearest == self.tscrunch)]
    
    def create_dm_list(self):
        DM_values = [pulsar['DM'] for pulsar in self.inj_report['pulsars']]
        DD_plan = self.create_DDplan()

        dm_trials_list = []
        tscrunch_list = []
        for i, dm_range in enumerate(DD_plan):
            n_trial = int(round((dm_range.high_dm-dm_range.low_dm)/dm_range.dm_step))
            if (i == 2) or (i == 0):
                n_trial, endpoint = (n_trial, False)
            else:
                n_trial, endpoint = (n_trial + 1, True) 


            dm_list_i = np.linspace(dm_range.low_dm, dm_range.high_dm, n_trial, endpoint=endpoint)
            tscrunch_i = np.ones_like(dm_list_i) * dm_range.tscrunch

            dm_trials_list.append(dm_list_i)
            tscrunch_list.append(tscrunch_i)
        
        dm_trials_arr = np.concatenate(dm_trials_list)
        tscrunch_arr = np.concatenate(tscrunch_list)
        dm_search_values = [self.find_n_nearest(dm_trials_arr, tscrunch_arr, DM_i, self.n_nearest) for DM_i in DM_values]
        dm_search_arr = np.unique(np.concatenate(dm_search_values))

        outfile = os.path.join(self.work_dir, f"dm_list_t{int(self.tscrunch)}.ascii")
        np.savetxt(outfile, dm_search_arr, fmt='%.3f')
        return outfile
    
    def run_cmd(self):
        chan_mask_file, birdie_list_file = self.generate_files()
        dm_list = self.create_dm_list()
        fft_size = int(round(self.processing_args['fft_length'] // self.tscrunch))
        ram_limit = self.processing_args['ram_limit']

        cmd = f"peasoup -k {chan_mask_file} -z {birdie_list_file} -i {self.fb} --dm_file {dm_list} " \
              f"--limit {self.processing_args['candidate_limit']} -n {self.processing_args['nharmonics']}  -m {self.processing_args['snr_threshold']} " \
              f"--acc_start {self.processing_args['start_accel']} --acc_end {self.processing_args['end_accel']} --fft_size {fft_size} " \
              f"-o {self.work_dir} --ram_limit_gb {ram_limit} --dedisp_gulp {self.gulp_size}"

        subprocess.run(cmd, shell=True)

    def peasoup_success(self, xml_name_old, xml_name_new, peasoup_dir):
        cand_finder = CandFinder(self.fb, self.inj_file, self.inj_report_path)
        cands_df = cand_finder.parse_xml_file(xml_name_old)
        self.n_harmonics = 2
        for n in range(self.n_harmonics):
            cands_data = cand_finder.filter_df(cands_df, snr_limit=5, pfact=n+1, adjust=0.05)
            cands_data['tscrunch'] = self.tscrunch
            cands_data.to_csv(f'{self.work_dir}/injected_xml_candidates_harm_{n+1}.csv')
        
        subprocess.run(f"rsync -Pav {xml_name_old} {peasoup_dir}/{xml_name_new}", shell=True)
        for n in range(self.n_harmonics):
            cand_name_old = f'injected_xml_candidates_harm_{n+1}.csv'
            cand_name_new = f'injected_xml_candidates_harm_{n+1}_tscrunch_{self.tscrunch}.csv'
            subprocess.run(f"rsync -Pav {self.work_dir}/{cand_name_old} {peasoup_dir}/{cand_name_new}", shell=True)
    
    def peasoup_failed(self, xml_name_new, xml_name, peasoup_dir):
        pointing_id = re.search(r'MMGPS_U_\d{4}', self.inj_report['injection_report']['fb']).group()
        inj_beam_name = re.search(r'[ci]fbf\d{5}', self.inj_report['injection_report']['fb']).group()
        xml_original = f'{self.data_dir}/{pointing_id}/XML_FILES/{inj_beam_name}/{xml_name}'
        subprocess.run(f"rsync -Pav {xml_original} {peasoup_dir}/{xml_name_new}", shell=True)

    def transfer_products(self):
        DD_plan = self.create_DDplan()
        xml_name = [f'overview_dm_{dm_range.low_dm:.6f}_{dm_range.high_dm:.6f}.xml' for dm_range in DD_plan][self.tscrunch_index]

        inj_ID = self.inj_report['injection_report']['ID']
        xml_name_old = f'{self.work_dir}/overview.xml'
        xml_name_new = f'{self.data_ID}_{inj_ID}_{xml_name}'
        peasoup_dir = f'{self.out_dir}/inj_{self.injection_number:06}/processing'
        
        if os.path.exists(xml_name_old):
            self.peasoup_success(xml_name_old, xml_name_new, peasoup_dir)
            process_dir = f'{self.out_dir}/inj_{self.injection_number:06}/processing'
            os.remove(f"{process_dir}/{self.fb}")
        else:
            self.peasoup_failed(xml_name_new, xml_name, peasoup_dir)

            
if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='peasoup for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--tscrunch_index', metavar='file', required=True, help='tscrucnh of filterbank file to search')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_file', metavar='file', required=True, help='JSON file with injection plan')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    parser.add_argument('--data_dir', metavar='dir', required=True, help='directory of xml files')


    args = parser.parse_args()

    peasoup_exec = PeasoupExec(args.tscrunch_index, args.search_args, args.injection_file, args.out_dir, args.data_dir, args.injection_number, -1)
    peasoup_exec.run_cmd()
    peasoup_exec.transfer_products()