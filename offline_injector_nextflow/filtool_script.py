import os
import sys
import json
import argparse
import subprocess
from collections import namedtuple


class FiltoolExec:
    def __init__(self, fb, search_args, output_dir, num_threads):
        self.fb = self.parse_fb(fb)
        self.num_threads = num_threads
        self.out = output_dir

        self.processing_args, self.ID = self.parse_search_args(search_args)
        self.define_defaults()
    
    @staticmethod
    def parse_fb(fb):
        if type(fb) == list:
            return ' '.join(fb)
        else:
            return fb

    def parse_search_args(self, search_flags):
        try:
            with open(search_flags, 'r') as file:
                flags = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {search_flags}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {search_flags} using JSON.')

        return flags['processing_args'], flags['processing_id']
    
    def define_defaults(self):
        self.pars = {
            'zapping_threshold': 4, 
            'segment_length': 2, 
            'baseline': '0 0', 
            'fillPatch': 'rand',
            'nbits': 8,
            'outmean': 128,
            'outstd': 6,
            'frequency_downsample': 1,
            'baseline_width': 0
        }

    def create_DDplan(self):
        DMRange = namedtuple("DMRange", ["low_dm", "high_dm", "dm_step", "tscrunch"])

        segments = []
        plan = self.processing_args['ddplan']
        for line in plan.splitlines():
            low_dm, high_dm, dm_step, tscrunch = list(map(float, line.split()[:4]))
            segments.append(DMRange(low_dm, high_dm, dm_step, tscrunch))

        return iter(sorted(segments, key=lambda x: x.tscrunch))
    
    def create_filplan(self, ddplan):
        filplan_file = os.path.join(self.out, "filplan.json")

        with open(filplan_file, 'w') as filplan:
            plans = []
            for dm_range in ddplan:
                plans.append({"time_downsample": dm_range.tscrunch,
                            "frequency_downsample": self.pars['frequency_downsample'],
                            "baseline_width": self.pars['baseline_width'],
                            "dataout_mean": self.pars['outmean'],
                            "dataout_std": self.pars['outstd'],
                            "dataout_nbits": self.pars['nbits'],
                            "rfi_flags": ""})
            json.dump(plans, filplan, indent=4)
        return filplan_file

    def run_cmd(self):
        ddplan = self.create_DDplan()
        filplan_file = self.create_filplan(ddplan)

        fscrunch = self.processing_args.get('fscrunch', 1) 
        rfi_flags = self.processing_args.get('rfi_flags', 'zdot')
        rootname = f"temp_merge_p_id_{self.ID}" 

        cmd = f"filtool -v -t {self.num_threads} --zapthre {self.pars['zapping_threshold']} --baseline {self.pars['baseline']} -l {self.pars['segment_length']} \
                 --filplan {filplan_file} --fd {fscrunch} --fillPatch {self.pars['fillPatch']} -z {rfi_flags} -o {self.out}/{rootname} -f {self.fb}"
        
        subprocess.run(cmd, shell=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='filtool for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--fb', metavar='file', required=True, nargs='+', help='injected filterbank file(s)')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    parser.add_argument('--n_threads', metavar='int', type=int, required=True, help='number of threads to use')
    args = parser.parse_args()

    fil_exec = FiltoolExec(args.fb, args.search_args, args.output, args.n_threads)
    fil_exec.run_cmd()

