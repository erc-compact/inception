import re
import sys
import json
import argparse
import subprocess

class FoldExec:
    def __init__(self, fb, search_args, injection_report, output, num_threads):
        self.out = output
        self.fb = self.parse_fb(fb)
        self.num_threads = num_threads
        args = self.parse_JSON(search_args)
        self.fold_args = args['fold_args']

        self.inj_report = self.parse_JSON(injection_report)
    
    @staticmethod
    def parse_fb(fb):
        if type(fb) == list:
            return ' '.join(fb)
        else:
            return fb

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
        
    def create_cand_file(self):
        pass
        
    def create_zap_sting(self):
        cmask = self.fold_args['channel_mask']
        cmask = cmask.strip()
        zap_string = ' '.join(['--rfi zap {} {}'.format(*i.split(':')) for i in cmask.split(',')])
        return zap_string
    
    def get_beam_tag(self):
        beam_name = re.search(r'[ci]fbf\d{5}', self.inj_report['injection']['fb']).group()
        if 'ifbf' in beam_name:
            beam_tag = "--incoherent"
        elif 'cfbf' in beam_name:
            beam_tag = f"-i {int(beam_name.strip('cfbf'))}"
        return beam_tag

    def run_cmd(self):
        TEMPLATE = '/home/psr/software/PulsarX/include/template/meerkat_fold.template'
        cand_file = self.create_cand_file()
        zap_string = self.create_zap_sting()
        beam_tag = self.get_beam_tag()

        nsubband = self.fold_args.get('nsubband', 64)
        fast_nbins = self.fold_args.get('fast_nbins', 64)
        slow_nbins = self.fold_args.get('slow_nbins', 128)

        cmd = f"psrfold_fil2 --dmboost 250 --plotx -v -t {self.num_threads} --candfile {cand_file} -n {nsubband} {beam_tag} " \
              f"-b {fast_nbins} --nbinplan 0.1 {slow_nbins} --template {TEMPLATE} --clfd 8 -L {self.fold_args['subint_length']} --fillPatch rand " \
              f"-f {self.fb} --rfi zdot {zap_string} --fd {self.fold_args['fscrunch']} --td {self.fold_args['tscrunch']}"

        subprocess.run(cmd, shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='candidate folder for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--fb', metavar='file', required=True, nargs='+', help='filterbank(s) to fold')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_report', metavar='file', required=True, help='JSON file with inject pulsar records')
    parser.add_argument('--output', metavar='dir', required=True, help='output directory')
    parser.add_argument('--n_threads', metavar='int', type=int, required=True, help='number of threads to use')
    args = parser.parse_args()

    fold_exec = FoldExec(args.fb, args.search_args, args.injection_report, args.output, args.n_threads)
    fold_exec.run_cmd()



"pred_file = generate_pulsarX_cand_file(
                        tmp_dir,
                        beam_name,
                        utc_start,
                        cand_mod_periods,
                        cand_dms,
                        cand_accs,
                        cand_snrs,
                        batch_start,
                        batch_stop)"