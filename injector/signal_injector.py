import os
import sys
import numpy as np 
from time import time
from pathlib import Path
from multiprocessing import Manager, Pool
from scipy.stats import truncnorm

from .io_tools import FilterbankReader, FilterbankWriter, print_exe
from .binary_model import BinaryModel
from .pulsar_model import PulsarModel
from .observation import Observation



class InjectSignal:
    def __init__(self, setup_manager, n_cpus, gulp_size_GB=0.01):
        manager = Manager()
        self.bits_flipped = manager.list([0]*n_cpus)  

        self.n_cpus = n_cpus
        self.gulp_size_GB = gulp_size_GB
        self.load_fb_stats = [setup_manager.fb.fb_mean, setup_manager.fb.fb_std]
        self.n_samples = setup_manager.fb.n_samples
        self.nchans = setup_manager.fb.nchans
        self.nbits = setup_manager.fb.nbits
        if (self.nchans * self.nbits) % 8:
            sys.exit(f'nchans x nbits ({self.nchans} x {self.nbits}) must be a multiple of 8.')
        self.signal_floor = 0.1 / (self.n_samples * self.nchans)
        self.compute_plan = self.create_parallel_plan()

        self.fb_path = setup_manager.fb.path
        self.seed = setup_manager.seed
        self.ephem = setup_manager.ephem
        self.out_path = setup_manager.output_path
        self.pulsars = setup_manager.pulsars
        self.parfile_paths = setup_manager.parfile_paths
        self.injected_path = self.out_path + '/' + Path(self.fb_path).stem + '_' + setup_manager.inj_ID
        self.part_path = self.injected_path + '.fil.part'

    def create_parallel_plan(self):
        block_size = int(self.gulp_size_GB/(self.nchans * self.nbits * 1.25e-10 * 2))
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
        for cpu_i in range(cpu):
            cpu_start += sum_file(self.compute_plan[cpu_i])

        return cpu_start
            
    def create_output(self):
        writer = FilterbankWriter(self.fb_path, self.part_path)
        data_bytes = os.path.getsize(self.fb_path) - writer.fb_reader.read_data_pos
        writer.write_file.truncate(writer.write_data_pos + data_bytes)
        writer.write_file.close()
        writer.fb_reader.read_file.close()

    def open_output_fb(self, cpu):
        filterbank_reader = FilterbankReader(self.fb_path, load_fb_stats=self.load_fb_stats)
        start = self.get_file_start(cpu)
        filterbank_reader.read_file.seek(filterbank_reader.read_data_pos + start*filterbank_reader.nchans*filterbank_reader.nbits//8)
        return FilterbankWriter(filterbank_reader, self.part_path, sample_offset=start)
    
    def construct_models(self, fb, cpu):
        pulsar_models = []
        for pulsar_data in self.pulsars:
            obs = Observation(fb, self.ephem, pulsar_data, generate=True)
            binary = BinaryModel(pulsar_data, generate=True)
            pulsar_model = PulsarModel(obs, binary, pulsar_data, generate=True)
            pulsar_models.append(pulsar_model)

        return pulsar_models

    def de_digitize(self, fb, values, rng):
        lower = (values - 0.5 - fb.fb_mean) / fb.fb_std
        upper = (values + 0.5 - fb.fb_mean) / fb.fb_std
        return truncnorm.ppf(rng.random(values.shape), lower, upper, loc=fb.fb_mean, scale=fb.fb_std)
    
    def inject_block(self, filterbank, cpu, block_start, block_size, models, rng):
        reader = filterbank.fb_reader
        block = reader.read_block(block_size)
        sample_start = block_start + self.get_file_start(cpu)
        
        pulsar_signal = np.zeros_like(block)
        for pulsar_model in models:
            pulsar_signal += pulsar_model.generate_signal(block_size, sample_start)

        channel_sigma = np.std(block, axis=0)
        pulsar_signal.T[channel_sigma==0] = 0

        active = np.abs(pulsar_signal) > self.signal_floor
        injected_block = block.copy()
        injected_block[active] = np.round(self.de_digitize(reader, block[active], rng) + pulsar_signal[active])
        filterbank.write_block(injected_block)

        self.bits_flipped[cpu] += int(np.count_nonzero(np.clip(injected_block, 0, 2**self.nbits-1) != block))


    def progress(self, cpu, N_blocks, block_i, t_stamp):
        if (block_i%10 == 0) and (block_i!=0):
            print_exe(f"CPU {cpu} processing {N_blocks} blocks: 10 blocks processed in {time()-t_stamp:.1f} s, {N_blocks-block_i} blocks remaining...")
            t_stamp = time()
        return t_stamp

    def inject_signal(self, cpu):
        fb = self.open_output_fb(cpu)
        models = self.construct_models(fb.fb_reader, cpu)
        rng = np.random.default_rng([self.seed, cpu])
        print_exe('Models constructed, starting injection...') if cpu == 0 else None
        (N_L_blocks, size_L_blocks), (_, size_S_blocks) = self.compute_plan[cpu]

        t_stamp = time()
        for block_i in range(N_L_blocks): 
            t_stamp = self.progress(cpu, N_L_blocks+int(size_S_blocks != 0), block_i, t_stamp)
            self.inject_block(fb, cpu, block_i*size_L_blocks, size_L_blocks, models, rng)

        if size_S_blocks != 0:
            self.inject_block(fb, cpu, N_L_blocks*size_L_blocks, size_S_blocks, models, rng)

        fb.fb_reader.read_file.close()
        fb.write_file.close()

    def parallel_inject(self):
        self.create_output()
        with Pool(self.n_cpus) as p:
            p.map(self.inject_signal, range(self.n_cpus))

    def finalise_output(self):
        os.replace(self.part_path, self.injected_path + '.fil')

        n_bits_flipped = np.sum(list(self.bits_flipped))
        print(f'bits flipped: {n_bits_flipped}/{self.n_samples*self.nchans} ({n_bits_flipped/(self.n_samples*self.nchans)*100:.3f}%)')