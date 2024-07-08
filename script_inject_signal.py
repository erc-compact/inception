import os
import sys
import getopt
import numpy as np
from time import time
from pathlib import Path
from multiprocessing import Pool
from scipy.stats import truncnorm
from inception.filterbank_tools import FilterbankIO
from inception.observatory import Observation
from inception.binary_model import PulsarBinaryModel
from inception.signal_generator import PulsarSignal


def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + output)

def progress(i, cpu, n_blocks, t_stamp):
    if (i%10 == 0) and (i!=0):
        print_exe(f"CPU {cpu} processing {n_blocks} blocks: 10 blocks processed in {time()-t_stamp:.1f} s, {n_blocks-i} blocks remaining...")
        t_stamp = time()
    return t_stamp

def de_digitize(fb_data, mean, sigma):

    def get_rvs(val):
        centre = (val-mean)/sigma
        deviation = 0.5/sigma
        d_plus = centre + deviation
        d_minus = centre - deviation
        return truncnorm(a=min(d_plus, d_minus), b=max(d_plus, d_minus), loc=mean, scale=sigma).rvs

    def de_digitizing(val):
        inds = np.where(fb_data == val)
        sampler = get_rvs(val)
        fb_data[inds] = sampler(size=len(inds[0]))
    
    for data in range(int(fb_data.min()), int(fb_data.max())+1):
        de_digitizing(data)

    return fb_data

def get_model(filterbank):
    obs = Observation(filterbank)
    binary = PulsarBinaryModel(pulsar_mass=pulsar_pars['Mass_p'], companion_mass=pulsar_pars['Mass_c'],
                               period=pulsar_pars['binary_period'], eccentricity=pulsar_pars['ecc'], 
                               inclination=pulsar_pars['inc'], LoAN=pulsar_pars['LoAN'], AoP=pulsar_pars['AoP'])
    pulsar_model = PulsarSignal(pulsar_pars['pulse_period'], pulsar_pars['duty_cycle'],
                                pulsar_pars['spectral_index'], pulsar_pars['SNR'], pulsar_pars['DM'],
                                obs, binary)
    return pulsar_model

def inject_signal(i, cpu, filesize, block_size, filterbank, pulsar_model):
    block = filterbank.read_block(block_size)
    pulsar_signal = pulsar_model.generate_signal(block_size, i*block_size + cpu*filesize)
    mean, sigma = filterbank.fb_mean, filterbank.fb_std
    analog_fb = de_digitize(block, mean, sigma)
    filterbank.write_block(np.round(analog_fb + pulsar_signal))

def create_files(cpu, filesize, n_blocks, block_size):
    injected_filepath = ad['--output'] + "/" + Path(ad['--filterbank']).stem + \
                        f"_{pulsar_pars['name']}"
    
    filterbank_sub = FilterbankIO(ad['--filterbank']) 
    filterbank_sub.new_filterbank(injected_filepath+f"_{cpu}.tmpfil")
    filterbank_sub.read_file.seek(filterbank_sub.read_data_pos + cpu*filesize*filterbank_sub.header['nchans'])
 
    pulsar_model = get_model(filterbank_sub)

    t_stamp = time()
    for i in range(n_blocks): 
        t_stamp = progress(i, cpu, n_blocks, t_stamp)
        inject_signal(i, cpu, filesize, block_size, filterbank_sub, pulsar_model)
    filterbank_sub.read_file.close()
    filterbank_sub.write_file.close()


def collect_files(n_files, n_blocks, block_size):
    injected_filepath = ad['--output'] + "/" + Path(ad['--filterbank']).stem + \
                        f"_{pulsar_pars['name']}"
    
    filterbank_main = FilterbankIO(ad['--filterbank']) 
    filterbank_main.new_filterbank(injected_filepath+".fil")
      
    for nf in range(n_files):
        print_exe(f'collecting file {nf+1}/{n_files}...')
        filterbank_sub = FilterbankIO(injected_filepath + f"_{nf}.tmpfil") 
        filterbank_sub.read_file.seek(filterbank_sub.read_data_pos)
        for _ in range(n_blocks):
            sub_block = filterbank_sub.read_block(block_size)
            filterbank_main.write_block(sub_block)
        filterbank_sub.read_file.close()
    filterbank_main.read_file.close()
    filterbank_main.write_file.close()


if __name__=='__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v", ["signal=", "filterbank=", "output="])
    except:
        print_exe(sys.exc_info()[0])
        sys.exit(2)

    ad = dict(opts) 
    pulsar_pars = dict(np.genfromtxt(ad['--signal'], dtype=None, encoding=None))
    keys = list(pulsar_pars.keys())
    for i in range(1,len(pulsar_pars)):
        key = keys[i]
        pulsar_pars[key] = np.float64(pulsar_pars[key])

    filterbank = FilterbankIO(ad['--filterbank'])   
    obs = Observation(filterbank)
    
    n_cpus = 2**5
    block_size = 2**11
    filesize = int(np.ceil(obs.n_samples/n_cpus))
    n_blocks = int(np.ceil(filesize/block_size))
    args = [(i, filesize, n_blocks, block_size) for i in range(n_cpus)]
    with Pool(n_cpus) as p:
        p.starmap(create_files, args)


    collect_files(n_cpus, n_blocks, block_size)