import os
import sys
import json
import argparse
import subprocess
import numpy as np
from pathlib import Path
from collections import namedtuple

from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

from inception.injector.io_tools import FilterbankReader


class PeasoupExec:
    def __init__(self, fb, search_args, injection_report, output_dir, ram_limit_gb, n_nearest):
        self.fb = fb
        self.out = output_dir
        self.ram_limit_gb = ram_limit_gb
        self.n_nearest = n_nearest

        args = self.parse_JSON(search_args)
        self.processing_args, self.ID  = args['processing_args'], args['processing_id']
        self.inj_report = self.parse_JSON(injection_report)

        self.gulp_size = self.get_gulp_size()
        self.tscrunch = self.get_tscrunch()
    
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
    
    def get_gulp_size(self):
        fscrunch = int(self.processing_args.get('fscrunch', 1))

        fb_reader = FilterbankReader(self.fb)
        default_gulp_size = int((2048.0 / (fb_reader.nchans / fscrunch)) * 1e6)
        del fb_reader
        return self.processing_args.get('gulp_size', default_gulp_size)
    
    def get_tscrunch(self):
        tscrunch_index = int(Path(self.fb).stem[-1]) - 1
        DD_plan = self.create_DDplan()
        
        return [dm_range.tscrunch for dm_range in DD_plan][tscrunch_index]
    
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
   
    def create_DDplan(self):
        DMRange = namedtuple("DMRange", ["low_dm", "high_dm", "dm_step", "tscrunch"])

        segments = []
        plan = self.processing_args['ddplan']
        for line in plan.splitlines():
            low_dm, high_dm, dm_step, tscrunch = list(map(float, line.split()[:4]))
            segments.append(DMRange(low_dm, high_dm, dm_step, tscrunch))

        return list(sorted(segments, key=lambda x: x.tscrunch))
    
    def generate_files(self):
        chan_mask_csv = self.processing_args['channel_mask']
        chan_mask_file = os.path.join(self.out, 'channel_mask.ascii')
        self.generate_chan_mask(chan_mask_csv, chan_mask_file)

        birdie_list_csv = self.processing_args['birdie_list']
        birdie_list_file = os.path.join(self.out, 'birdie_list.ascii')
        self.generate_birdie_list(birdie_list_csv, birdie_list_file)
        return chan_mask_file, birdie_list_file
    
    @staticmethod
    def find_n_nearest(dm_arr, tscrunch_arr, target, n):
        if n == -1:
            return dm_arr, tscrunch_arr
        else:
            nearest_indices = np.argsort(np.abs(dm_arr - target))[:n]
            dm_values = dm_arr[nearest_indices]
            tscrunch_values = tscrunch_arr[nearest_indices]
            return dm_values, tscrunch_values
    
    def create_dm_list(self):
        DM_values = [pulsar['DM'] for pulsar in self.inj_report['pulsars']]
        DD_plan = self.create_DDplan()

        dm_trials_list = []
        tscrunch_list = []
        for i, dm_range in enumerate(DD_plan):
            n_trial = int((dm_range.high_dm-dm_range.low_dm)/dm_range.dm_step)
            n_trial, endpoint = (n_trial + 1, True) if (i == len(DD_plan)-1) else (n_trial, False)

            dm_list_i = np.linspace(dm_range.low_dm, dm_range.high_dm, n_trial, endpoint=endpoint)
            tscrunch_i = np.ones_like(dm_list_i) * dm_range.tscrunch

            dm_trials_list.append(dm_list_i)
            tscrunch_list.append(tscrunch_i)
        
        dm_trials_arr = np.concatenate(dm_trials_list)
        tscrunch_arr = np.concatenate(tscrunch_list)

        dm_search_values = []
        for DM_i in DM_values:
            dm_nearest, tscrunch_nearest = self.find_n_nearest(dm_trials_arr, tscrunch_arr, DM_i, self.n_nearest)
            dm_search_values.append(dm_nearest[np.where(tscrunch_nearest == self.tscrunch)])

        dm_search_arr = np.concatenate(dm_search_values)
        dm_search_arr = np.unique(dm_search_arr)

        outfile = os.path.join(self.out, f"dm_list_t{int(self.tscrunch)}.ascii")
        np.savetxt(outfile, dm_search_arr, fmt='%.3f')
        return outfile
    
    def rename_outputs(self):
        DD_plan = self.create_DDplan()
        tscrunch_index = int(Path(self.fb).stem[-1]) - 1
        xml_name = [f'overview_dm_{dm_range.low_dm:.6f}_{dm_range.high_dm:.6f}.xml' for dm_range in DD_plan][tscrunch_index]

        inj_ID = self.inj_report['injection']['ID']
        xml_name_new = f'{self.out}/{inj_ID}_{self.ID}_{xml_name}'

        os.rename(f'{self.out}/overview.xml', xml_name_new)
    
    def run_cmd(self):
        chan_mask_file, birdie_list_file = self.generate_files()
        dm_list = self.create_dm_list()

        cmd = f"peasoup -k {chan_mask_file} -z {birdie_list_file} -i {self.fb} --dm_file {dm_list} " \
              f"--limit {self.processing_args['candidate_limit']} -n {self.processing_args['nharmonics']}  -m {self.processing_args['snr_threshold']} " \
              f"--acc_start {self.processing_args['start_accel']} --acc_end {self.processing_args['end_accel']} --fft_size {self.processing_args['fft_length']}" \
              f" -o {self.out} --ram_limit_gb {self.ram_limit_gb} --dedisp_gulp {self.gulp_size}"
        
        subprocess.run(cmd, shell=True)

        self.rename_outputs()


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='peasoup for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--fb', metavar='file', required=True, help='Injected filterbank file')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_report', metavar='file', required=True, help='JSON file with inject pulsar records')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    parser.add_argument('--ram_limit', metavar='GB', required=True, help='limit the ram used by a single peasoup search')
    parser.add_argument('--n_nearest', metavar='int', type=int, required=True, help='number of DM trials to search around injected pulsar DM')
    args = parser.parse_args()

    peasoup_exec = PeasoupExec(args.fb, args.search_args, args.injection_report, args.output, args.ram_limit, args.n_nearest)
    peasoup_exec.run_cmd()
    