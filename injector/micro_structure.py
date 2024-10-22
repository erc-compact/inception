import sys
import numpy as np 
import astropy.units as u


class MicroStructure: 
    def __init__(self, phase_abs, freq, scale, pulse_period, profile, seed):
        self.scale = pulse_period/(scale*u.us.to(u.s))
        self.seed = seed

        self.phase_abs = phase_abs
        self.intrinsic_profile = profile(phase_abs % 1, freq)
        self.noise = np.zeros_like(phase_abs)
        self.n_samples, self.nchans = phase_abs.shape

        self.pulse_numbers = np.floor(phase_abs).astype(int)
        self.pulse_range = np.arange(np.min(self.pulse_numbers), np.max(self.pulse_numbers)+1)

        self.pulse_counts_block = {}
        self.pulse_index = {}

        self.max_pulse_length = self.pre_process()

        self.profile = self.pulse_profile()

    def pulse_sample_counter(self, pulse_number):
        pulse_index = np.where(self.pulse_numbers == pulse_number)
        pulse_counts = np.bincount(pulse_index[1], minlength=self.nchans)
        return pulse_counts, pulse_index

    @staticmethod
    def smoothstep_S2(x):
        return 6 * x**5 - 15 * x**4 + 10 * x**3
    
    @staticmethod
    def get_pulse_rng(pulse_num, seed=0):
        pulse_offset = 10**9 if np.sign(pulse_num) == -1 else 0
        return np.random.default_rng(seed+pulse_num+pulse_offset)
    
    def pre_process(self):
        
        for pulse_n in self.pulse_range:
            pcb, pix = self.pulse_sample_counter(pulse_n)
            self.pulse_counts_block[pulse_n] = pcb
            self.pulse_index[pulse_n] = pix
        
        max_pulse_length = np.max([np.max(pcb) for pcb in self.pulse_counts_block.values()])
        return max_pulse_length

    def perlin_noise(self, pulse_length, pulse_num):
        rng = self.get_pulse_rng(pulse_num, self.seed)
        scale = abs(rng.normal(self.scale, self.scale*0.05))
        gradients = rng.uniform(-1, 1, int(scale) + 1)
        gradients /= np.linalg.norm(gradients, axis=0)  

        grid_points = np.arange(pulse_length)
        norm_len = pulse_length / scale
        x0 = grid_points / norm_len
        rel_x = (grid_points % norm_len) / norm_len
        dot0 = gradients[(x0 % scale).astype(int)] * rel_x
        dot1 = gradients[((x0+1) % scale).astype(int)] * (rel_x - 1)

        u = self.smoothstep_S2(rel_x)
        noise = (1 - u) * dot0 + u * dot1
        return np.abs(noise)
    
    @staticmethod
    def get_pad_chunks(pad_arr):
        nonzero = pad_arr != 0
        changes = np.diff(nonzero.astype(int))

        chunk_starts = np.where(changes == 1)[0] + 1
        chunk_ends = np.where(changes == -1)[0] + 1

        if nonzero[0]:
            chunk_starts = np.r_[0, chunk_starts]

        if nonzero[-1]:
            chunk_ends = np.r_[chunk_ends, len(pad_arr)]

        return chunk_starts, chunk_ends
    
    @staticmethod
    def get_padding(chan, pad, chunk_starts, chunk_ends):
        chan_pad = pad[chan]
        if (chan_pad == 0) or (len(chunk_starts) == 0):
            s_pad, e_pad = 0, None

        elif len(chunk_starts) == 2:
            if chunk_starts[0] <= chan <  chunk_ends[0]:
                s_pad, e_pad = chan_pad, None
            elif chunk_starts[1] <= chan <  chunk_ends[1]:
                s_pad, e_pad = 0, -chan_pad

        elif len(chunk_starts) == 1:
            if pad[chunk_starts[0]] > pad[chunk_ends[0]-1]:
                s_pad, e_pad = chan_pad, None
            else:
                s_pad, e_pad = 0, -chan_pad
        else:
            sys.exit("I was aware of this edge-case but didn't think it was possible. Please let me know what inputs casued this error.")
        
        return s_pad, e_pad
    
    def create_microstructure(self, pulse_num):
        pulse_counts_block = self.pulse_counts_block[pulse_num]
        pulse_index = self.pulse_index[pulse_num]
   
        pulse_counts = pulse_counts_block.copy()
        max_pulse_length = self.max_pulse_length 
        pulse_counts[(pulse_counts_block < max_pulse_length-2) & (pulse_counts_block > 0)] = max_pulse_length # fix edge case
        
        pad = pulse_counts - pulse_counts_block
        chunk_starts, chunk_ends = self.get_pad_chunks(pad)

        noise_unique = {max_pulse_length: self.perlin_noise(max_pulse_length, pulse_num), #TODO: do this recursively
                        max_pulse_length-1: self.perlin_noise(max_pulse_length-1, pulse_num),
                        max_pulse_length-2: self.perlin_noise(max_pulse_length-2, pulse_num)}


        for chan in range(self.nchans):
            pulse_intrinsic_length = pulse_counts[chan]
            if pulse_intrinsic_length != 0:
                s_pad, e_pad = self.get_padding(chan, pad, chunk_starts, chunk_ends)

                micro_structure_values = noise_unique[pulse_intrinsic_length]                

                noise_values = micro_structure_values[s_pad: e_pad]
                channel_index = pulse_index[0][np.where(pulse_index[1] == chan)[0]]

                intrinsic_profile = self.intrinsic_profile[channel_index, chan]
                profile_sum = np.sum(intrinsic_profile*noise_values)
                if profile_sum != 0:
                    norm = np.sum(intrinsic_profile) / profile_sum
                else:
                    norm = 0

                self.noise[channel_index, chan] = noise_values * norm

    def pulse_profile(self):
        for pulse_n in self.pulse_range:  
            self.create_microstructure(pulse_n)

        return self.intrinsic_profile * self.noise