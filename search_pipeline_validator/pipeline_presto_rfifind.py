import os
import glob
import argparse
import subprocess
from pathlib import Path

import pipeline_tools as inj_tools


class RFIPrestoProcess:
    def __init__(self, processing_args, out_dir, work_dir, injection_number):

        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number

    def setup(self):
        self.get_injection_report()
        self.transfer_data()

    def get_injection_report(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        report_path = glob.glob(f'{results_dir}/report_*.json')[0]
        self.injection_report = inj_tools.parse_JSON(report_path)
        self.inj_id = self.injection_report['injection_report']['ID']

    def transfer_data(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'

        data = glob.glob(f"{results_dir}/*_{self.inj_id}.fil")[0]
        inj_tools.rsync(data, self.work_dir)

        self.data = f'{self.work_dir}/{Path(data).name}'

    def run(self, threads):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing'
        os.makedirs(presto_out_dir, exist_ok=True)

        rfi_args = self.processing_args['presto_rfi_args']

        cmd = f"rfifind -ncpus {threads} -time {rfi_args['rfifind']['time']} -o {presto_out_dir}/{self.inj_id} {self.data}"
        subprocess.run(cmd, shell=True)

        out_file=f"{self.work_dir}/{self.inj_id}_topo_DM0.00"
        cmd=f"prepdata -ncpus {threads} -nobary -o {out_file} -dm 0.0 -mask {presto_out_dir}/{self.inj_id}_rfifind.mask {self.data}"
        subprocess.run(cmd, shell=True)

        cmd=f"realfft {out_file}.dat"
        subprocess.run(cmd, shell=True)

        cmd=f"accelsearch -ncpus {threads} -numharm {rfi_args['birdies']['numharm']} -zmax 0 {out_file}.fft"
        subprocess.run(cmd, shell=True)

        subprocess.run("ls", shell=True)

    def make_birdies(self):
        pass


    def transfer_products(self):
        results_dir = f'{self.out_dir}/inj_{self.injection_number:06}'
        presto_out_dir = f'{results_dir}/processing'
        out_file=f"{self.work_dir}/{self.inj_id}_topo_DM0.00"

        inj_tools.rsync(f'{self.work_dir}/*.txt', presto_out_dir)
        inj_tools.rsync(f'{self.work_dir}/*.inf', presto_out_dir)
        inj_tools.rsync(f'{self.work_dir}/*.cand', presto_out_dir)
        inj_tools.rsync(f'{self.work_dir}/*ACCEL_0*', presto_out_dir)

        if self.processing_args['presto_rfi_args']['save_dat']:
            inj_tools.rsync(f'{out_file}.dat', presto_out_dir)

        if self.processing_args['presto_rfi_args']['save_fft']:
            inj_tools.rsync(f'{out_file}.dat', presto_out_dir)
            

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='presto rfi cleaning for search pipeline validator',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--injection_number', metavar='int', required=True, type=int, help='injection process number')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with search parameters')

    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    parser.add_argument('--threads', metavar='int', type=int, required=True, help='number of cpus to use')
    args = parser.parse_args()

    rfi_exec = RFIPrestoProcess(args.processing_args, args.out_dir, args.work_dir, args.injection_number)
    rfi_exec.setup()
    rfi_exec.run(args.threads)
    rfi_exec.make_birdies()
    rfi_exec.transfer_products()