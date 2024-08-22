import os
import sys
import numpy as np 
from time import time
from pathlib import Path
from multiprocessing import Pool
from scipy.stats import truncnorm

from injector.io_tools import FilterbankIO, print_exe
from injector.binary_model import PulsarBinaryModel
from injector.signal_generator import PulsarSignal
from injector.observatory import Observation


class InjectSignal:
    def __init__(self, setup_manager, n_cpus):
        self.n_cpus = int(n_cpus)
        self.n_samples = setup_manager.pulsar_models[0].obs.n_samples
        self.compute_plan = self.create_parallel_plan()

        self.fb_path = setup_manager.fb.path
        self.ephem = setup_manager.ephem
        self.out_path = setup_manager.output_path
        self.pulsars = setup_manager.pulsars
        self.parfile_paths = setup_manager.parfile_paths
        self.injected_path = self.out_path + '/' + Path(self.fb_path).stem + '_' + setup_manager.pulsar_models[0].ID

    def create_parallel_plan(self):
        block_size = 2**11 # optimise?
        filesize_per_cpu, filesize_remainder = divmod(self.n_samples, self.n_cpus)
        large_block_size, large_remainder = divmod(filesize_per_cpu+1, block_size)
        small_block_size, small_remainder = divmod(filesize_per_cpu, block_size)

        compute_plan = dict([(i, []) for i in range(self.n_cpus)])
        for cpu in range(self.n_cpus):
            if cpu < filesize_remainder:
                compute_plan[cpu].append((large_block_size, block_size))
                compute_plan[cpu].append((1, large_remainder))
            else:
                compute_plan[cpu].append((small_block_size, block_size))
                compute_plan[cpu].append((1, small_remainder))

        return compute_plan

    def get_file_start(self, cpu):
        def sum_file(cpu_info):
            (N_L_block, size_L_block), (N_S_block, size_S_block) = cpu_info
            return N_L_block*size_L_block + N_S_block*size_S_block

        cpu_start = 0
        for cpu in range(cpu):
            cpu_start += sum_file(self.compute_plan[cpu])

        return cpu_start
            
    def open_tmp_fb(self, cpu):        
        filterbank_sub = FilterbankIO(self.fb_path) 
        filterbank_sub.read_file.seek(filterbank_sub.read_data_pos + self.get_file_start(cpu)*filterbank_sub.header['nchans'])

        filterbank_sub.new_filterbank(self.injected_path + f"_{cpu}.tmpfil")
        return filterbank_sub
    
    def construct_models(self, fb):
        pulsar_models = []
        for i, pulsar_data in enumerate(self.pulsars):
            obs = Observation(fb, self.ephem, pulsar_data, validate=False)
            binary = PulsarBinaryModel(pulsar_data, validate=False)
            pulsar_model = PulsarSignal(obs, binary, pulsar_data, validate=False, par_file=self.parfile_paths[i])
            pulsar_models.append(pulsar_model)

        return pulsar_models

    @staticmethod
    def de_digitize(fb, data_block):

        def get_rvs(val):
            centre = (val-fb.fb_mean)/fb.fb_std
            deviation = 0.5/fb.fb_std
            d_plus = centre + deviation
            d_minus = centre - deviation
            return truncnorm(a=min(d_plus, d_minus), b=max(d_plus, d_minus), loc=fb.fb_mean, scale=fb.fb_std).rvs

        def de_digitizing(val):
            inds = np.where(data_block == val)
            sampler = get_rvs(val)
            data_block[inds] = sampler(size=len(inds[0]))
        
        for data in range(int(data_block.min()), int(data_block.max())+1):
            de_digitizing(data)

        return data_block
    
    def inject_block(self, fb, cpu, block_start, block_size, models):
        block = fb.read_block(block_size)
        sample_start = block_start + self.get_file_start(cpu)
        
        pulsar_signal = np.zeros_like(block)
        for pulsar_model in models:
            pulsar_signal += pulsar_model.generate_signal(block_size, sample_start)

        analog_block = self.de_digitize(fb, block)
        fb.write_block(np.round(analog_block + pulsar_signal))

    def progress(self, cpu, N_blocks, block_i, t_stamp):
        if (block_i%10 == 0) and (block_i!=0):
            print_exe(f"CPU {cpu} processing {N_blocks} blocks: 10 blocks processed in {time()-t_stamp:.1f} s, {N_blocks-block_i} blocks remaining...")
            t_stamp = time()
        return t_stamp

    def inject_signal(self, cpu):
        fb = self.open_tmp_fb(cpu)
        models = self.construct_models(fb)
        (N_L_blocks, size_L_blocks), (_, size_S_blocks) = self.compute_plan[cpu]

        t_stamp = time()
        for block_i in range(N_L_blocks): 
            t_stamp = self.progress(cpu, N_L_blocks+int(size_S_blocks != 0), block_i, t_stamp)
            self.inject_block(fb, cpu, block_i*size_L_blocks, size_L_blocks, models)

        if size_S_blocks != 0:
            self.inject_block(fb, cpu, N_L_blocks*size_L_blocks, size_S_blocks, models)

        fb.read_file.close()
        fb.write_file.close()

    def parallel_inject(self):
        args = list(range(self.n_cpus))

        with Pool(self.n_cpus) as p:
            p.map(self.inject_signal, args)

    def combine_files(self):
        filterbank_main = FilterbankIO(self.fb_path) 
        filterbank_main.new_filterbank(self.injected_path + ".fil")
        
        for cpu in range(self.n_cpus):
            print_exe(f'combining file {cpu+1}/{self.n_cpus}...')
            filterbank_sub = FilterbankIO(self.injected_path + f"_{cpu}.tmpfil") 
            filterbank_sub.read_file.seek(filterbank_sub.read_data_pos)

            (N_L_blocks, size_L_blocks), (_, size_S_blocks) = self.compute_plan[cpu]

            for _ in range(N_L_blocks):
                L_sub_block = filterbank_sub.read_block(size_L_blocks)
                filterbank_main.write_block(L_sub_block)

            S_sub_block = filterbank_sub.read_block(size_S_blocks)
            filterbank_main.write_block(S_sub_block)
            filterbank_sub.read_file.close()
            os.remove(self.injected_path + f"_{cpu}.tmpfil")

        filterbank_main.read_file.close()
        filterbank_main.write_file.close()

    

