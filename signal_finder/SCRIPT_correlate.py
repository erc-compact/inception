import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

import time
import getopt
import itertools
import numpy as np
from operator import itemgetter
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.signal import correlate, correlation_lags

from injector.observatory import Observation
from injector.io_tools import read_datfile, FilterbankIO, print_exe


class SignalCorrelate:
    def __init__(self, output, fb_name, ncpu=1, DM_min=0, DM_max=500, DM_resolution=2):
        self.output = output
        self.fb_name = fb_name
        self.ncpu = int(ncpu)

        self.DM_min = np.float64(DM_min)
        self.DM_max = np.float64(DM_max)
        self.DM_res = int(DM_resolution)

        self.DM_bins = np.round(np.arange(self.DM_min, self.DM_max, 10**-self.DM_res, dtype=np.float64), self.DM_res)
        self.DM_bins_dict = dict(zip(self.DM_bins, np.arange(len(self.DM_bins))))

        fb = FilterbankIO(fb_name) 
        fb.split_fb_channels(output)
        del fb

        self.stacked_power_cpu = {}

    def get_compute_data(self):
        fb = FilterbankIO(self.fb_name) 
        obs = Observation(fb, 'BUILTIN', {'DM':0, 'ID':'None'}, validate=True)
        channel_pairs = list(itertools.combinations(range(obs.n_chan), 2))
        channel_files = [read_datfile(self.get_path(chan), nbits=64) for chan in range(obs.n_chan)]
        return obs, channel_pairs, channel_files

    def get_path(self, channel):
        return self.output+f'/_{channel}.dat'

    @staticmethod
    def cross_corr(chan1, chan2, channel_files): 
        data1 = channel_files[chan1]
        data2 = channel_files[chan2]
        return correlate((data1-np.mean(data1)), (data2-np.mean(data2)), mode='same', method='fft')
    
    @staticmethod
    def dt2dm(chan1, chan2, obs):
        dt = correlation_lags(obs.n_samples, obs.n_samples, 'same')*obs.dt
        f1, f2 = obs.freq_arr[[chan1, chan2]]
        freq_div = (1/f2**2 - 1/f1**2)
        DM = -dt / (obs.DM_const * freq_div)
        return DM
    
    def get_correlation(self, channel_pair, obs, channel_files):
        power_stack = np.zeros_like(self.DM_bins)
  
        power_corr = self.cross_corr(*channel_pair, channel_files)
        DM_corr = self.dt2dm(*channel_pair, obs)

        DM_corr_rounded = np.round(DM_corr, self.DM_res)
        DM_range = (DM_corr_rounded > self.DM_min) & (DM_corr_rounded < self.DM_max)

        DM_corr_clipped = DM_corr_rounded[DM_range]

        bins = itemgetter(*DM_corr_clipped)(self.DM_bins_dict)
        power_stack[[*bins]] += power_corr[DM_range]
        return power_stack
    
    def stack_correlations(self, cpu, pair_inds):
        obs, channel_pairs, channel_files = self.get_compute_data()
        power_stack = np.zeros_like(self.DM_bins)

        compute_pairs = channel_pairs[slice(*pair_inds)]
        t_stamp = time.time()
        for i, pair in enumerate(compute_pairs):
            if (i % 50 == 0) and (i != 0):
                print_exe(f'cpu {cpu}: 50 cross-correlation pairs completed in {time.time()-t_stamp:.2f}s, {len(compute_pairs)-i} pairs remaining')
                t_stamp = time.time()

            power_stack += self.get_correlation(pair, obs, channel_files)

        self.stacked_power_cpu[cpu] = power_stack

    def compute_power_stack(self):
        _, channel_pairs, _ = self.get_compute_data()

        N_pairs, remain = divmod(len(channel_pairs), self.ncpu)
        start_pairs, end_pairs = [], []
        for cpu in range(self.ncpu):
            add = 1 if cpu < remain else 0 
            start_pairs.append(cpu*(N_pairs+add))
            end_pairs.append((cpu+1)*(N_pairs+add))

        pair_inds = zip(start_pairs, end_pairs)
        args = list(zip(np.arange(self.ncpu), pair_inds))
        with Pool(self.ncpu) as p:
            p.starmap(self.stack_correlations, args)

    def clean_directory(self):
        obs, _, _ = self.get_compute_data()
        for chan in range(obs.n_chan):
            os.remove(self.get_path(chan))

    def save_compute(self, save_plot=False, save_file=False):
        stacked_power = np.sum(np.stack(list(self.stacked_power_cpu.values())), axis=0)

        if save_file:
            np.save('stacked_DM', np.stack([self.DM_bins, stacked_power]))
        if save_plot:
            fig, ax = plt.subplots(figsize=(15, 4))
            ax.plot(self.DM_bins, stacked_power)
            ax.set_ylabel('stacked power')
            ax.set_xlabel('DM')

            fig.savefig(self.output+'/stacked_DM.png', bbox_inches="tight")
        
        return self.DM_bins, stacked_power
    

if __name__=='__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["filterbank=", "output=", 'ncpu=', 
                                                       "DM_min=", "DM_max=", "DM_res="])
    except Exception as err:
        sys.exit(err)
    ad = dict(opts) 

    computer = SignalCorrelate(ad['--output'], ad['--filterbank'], 
                               ncpu=ad.get('--ncpu', 1), 
                               DM_min=ad.get('--DM_min', 0), 
                               DM_max=ad.get('--DM_max', 500), 
                               DM_resolution=ad.get('--DM_res', 2))
    
    computer.compute_power_stack()
    computer.clean_directory()
    _ = computer.save_compute(save_plot=True, save_file=True)