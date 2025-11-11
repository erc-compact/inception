import os
import sys
import struct
import numpy as np 


class FilterbankReader:     
    def __init__(self, filterbank, gulp_size_GB=0.01, stats_samples=1e6, load_fb_stats=[]):
        self.inttypes = ['machine_id', 'telescope_id', 'data_type', 'nchans','nbits', 'nifs', 'scan_number', 
                         'barycentric','pulsarcentric', 'nbeams', 'ibeam', 'nsamples']
        self.strtypes =  ['source_name','rawdatafile']
        self.dbltypes = ['tstart','tsamp','fch1','foff','refdm','az_start','za_start','src_raj','src_dej']
        self.chrtypes = ['signed']

        self.header={}
        self.path = filterbank
        self.read_file = None
        self.read_header_pos = 0
        self.read_data_pos = 0
        self.output_dtype = np.float64
        self.gulp_size_GB = gulp_size_GB

        self.read_header(filterbank)
        self.dt = self.header['tsamp']
        self.nchans = int(self.header['nchans'])
        self.nbits = int(self.header['nbits'])
        self.bandwidth = abs(self.header['foff']) * self.header['nchans']
        self.ftop = self.header['fch1'] - 0.5 * self.header['foff']
        self.fbottom = self.ftop + self.header['foff'] * self.header['nchans']
        self.center = self.ftop + 0.5 * self.header['foff'] * self.header['nchans']

        self.n_samples = self.get_n_samples() 
        self.obs_len = self.n_samples * self.dt
        if load_fb_stats:
            self.fb_mean, self.fb_std = load_fb_stats
        else:
            if stats_samples:
                print_exe(f'calculating filterbank statistics using {int(stats_samples)}...')
                self.fb_mean, self.fb_std = self.get_FB_stats(stats_samples)
                print_exe(f'mean: {self.fb_mean}, std: {self.fb_std}')
            else:
                self.fb_mean, self.fb_std = 128.0, 6.0

    def read_string(self):
        nchar = np.fromfile(self.read_file, dtype=np.int32, count=1)[0]
        if nchar == 0:
            return ""
        if nchar < 1 or nchar > 80:
            sys.exit(f"Cannot parse filterbank header at {self.path}.")
        byte_data = self.read_file.read(nchar)
        string_data = byte_data.decode("UTF-8")
        return string_data
    
    def read_key_data(self, key):
        if key in self.inttypes:
            return np.fromfile(self.read_file, dtype=np.int32, count=1)[0]
        elif key in self.strtypes:
            return self.read_string()
        elif key in self.dbltypes:
            return np.fromfile(self.read_file, dtype=np.float64, count=1)[0]
        elif key in self.chrtypes:
            return np.fromfile(self.read_file,dtype=np.int8,count=1)[0]
        else:
            sys.exit(f"Cannot parse filterbank header, key '{key}' not understood")

    def read_header(self, filename):
        self.read_file = open(filename, 'rb')
        self.read_header_pos=self.read_file.tell()

        key = self.read_string()
        if key=="HEADER_START":
            key = self.read_string()
            while key != "HEADER_END":
                self.header[key] = self.read_key_data(key)
                key = self.read_string()
        else:
            sys.exit("Cannot parse filterbank header, HEADER_START was not found")
        self.read_data_pos=self.read_file.tell()

    # def get_FB_stats_gulp(self, nsamples):
    #     current_loc = self.read_file.tell()
    #     data = self.read_block(nsamples)
    #     self.read_file.seek(current_loc, 0)
    #     return np.median(np.mean(data, axis=0)), np.median(np.std(data, axis=0))
    
    @staticmethod
    def Welford_alg(count, mean, M2, new_data):
        n = new_data.shape[0]
        count += n

        delta = new_data - mean  
        new_mean = mean + np.sum(delta, axis=0) / count

        delta2 = new_data - new_mean
        M2 += np.sum(delta * delta2, axis=0)

        return count, new_mean, M2

    def get_FB_stats(self, stat_samples):
        count = 0
        mean = np.zeros(self.nchans)
        M2 = np.zeros(self.nchans)

        self.read_file.seek(self.read_data_pos, 0)
        n_chunks =  int(np.int64(self.n_samples) * self.nchans * self.nbits * 1.25e-10/self.gulp_size_GB)

        n_samples = int(stat_samples) if (stat_samples != 0) and (stat_samples < self.n_samples) else self.n_samples
        full_block_size, remainder = divmod(abs(n_samples), n_chunks)
        for chunk in range(n_chunks):
            if chunk < remainder:
                block_size = full_block_size + 1
            else:
                block_size = full_block_size

            self.read_file.seek(self.n_samples // n_chunks * self.nchans * chunk, 0)
            data = self.read_block(block_size)
            count, mean, M2 = self.Welford_alg(count, mean, M2, data)

        std = np.sqrt(M2 / (count - 1)) 
        global_mean = np.median(mean[std != 0])
        global_std = np.median(std[std != 0])
        self.read_file.seek(self.read_data_pos, 0)
        return global_mean, global_std
    
    def get_n_samples(self):
        self.read_file.seek(0,2)
        n_samples = int((self.read_file.tell() - self.read_data_pos)/self.nchans * (8/self.nbits))
        self.read_file.seek(self.read_data_pos)
        return n_samples

    def read_block(self, nsamples):
        nbytes = nsamples * self.nchans

        if self.nbits == 16:
            out = np.fromfile(self.read_file,dtype=np.uint16,count=nbytes)
        elif self.nbits == 8:
            out = np.fromfile(self.read_file,dtype=np.uint8,count=nbytes)
        elif self.nbits == 4:
            raw = np.fromfile(self.read_file, dtype=np.uint8, count=nbytes // 2)
            out = np.empty(nbytes, dtype=np.uint8)
            out[0::2] = raw & 0x0F 
            out[1::2] = (raw >> 4) & 0x0F 
        elif self.nbits == 2:
            raw = np.fromfile(self.read_file, dtype=np.uint8, count=nbytes // 4)
            out = np.empty(nbytes, dtype=np.uint8)
            out[0::4] = raw & 0x03         
            out[1::4] = (raw >> 2) & 0x03            
            out[2::4] = (raw >> 4) & 0x03             
            out[3::4] = (raw >> 6) & 0x03         
        elif self.nbits == 1:
            raw = np.fromfile(self.read_file, dtype=np.uint8, count=nbytes // 8)
            raw_bits = (raw[:, np.newaxis] >> np.arange(8)) & 1
            out = raw_bits.flatten()
        else:
            sys.exit(f"{self.nbits} bit filterbank is not supported.")
        
        return out.reshape(-1, self.nchans).astype(self.output_dtype)
    
    def split_fb_channels(self, output_path, ext='', nbits=32):
        block_size = 2**11
        n_blocks, remainder = divmod(self.n_samples, block_size)
        open_files = [open(output_path+f'/{ext}_{chan}.dat', 'wb') for chan in range(self.nchans)]
        data_type = {64: np.float64, 32: np.float32, 8: np.uint8}
        def split_block(block):
            for chan in range(self.nchans):
                raw = block.T[chan].astype(data_type[nbits])
                open_files[chan].write(raw.tobytes())

        for i in range(n_blocks):
            print_exe(f'{i}/{n_blocks} blocks processed') if (i%100 == 0) else None
            block = self.read_block(block_size)
            split_block(block)
        
        block_remain = self.read_block(remainder)
        split_block(block_remain)


class FilterbankWriter: 
    def __init__(self, read_filterbank, write_filterbank_name, gulp_size_GB=0.01, stats_samples=0, load_fb_stats=[]):
        self.fb_reader = FilterbankReader(read_filterbank, gulp_size_GB, stats_samples, load_fb_stats) if type(read_filterbank) == str else read_filterbank

        self.nbits = self.fb_reader.nbits
        self.write_file = None
        self.write_header_pos = 0 
        self.write_data_pos = 0

        self.create_filterbank(write_filterbank_name)
    
    def create_filterbank(self, filename):
        def write_string(file, string):
            l=len(string)
            byte_data= struct.pack(f"i{l}s", l, string.encode("UTF-8"))
            file.write(byte_data)
        
        self.write_file = open(filename, "wb")
        self.write_header_pos=self.write_file.tell()

        write_string(self.write_file,"HEADER_START")
        for key in self.fb_reader.header:
            write_string(self.write_file,key)
            if key in self.fb_reader.inttypes:
                self.write_file.write(struct.pack("i",self.fb_reader.header[key]))
            elif key in self.fb_reader.strtypes:
                write_string(self.write_file,self.fb_reader.header[key])
            elif key in self.fb_reader.dbltypes:
                self.write_file.write(struct.pack("d",self.fb_reader.header[key]))
            elif key in self.fb_reader.chrtypes:
                self.write_file.write(struct.pack("b",self.fb_reader.header[key]))
            else:
                raise Exception(f"Cannot understand filterbank header for writing, key '{key}' not understood")
        write_string(self.write_file,"HEADER_END")
        self.write_data_pos=self.write_file.tell()

    def write_block(self, block):
        block = np.clip(block, 0, 2**self.nbits-1)
        # block = np.mod(block, 2**self.nbits-1) overflows values, too slow

        if self.nbits == 16:
            out = block.flatten(order='C').astype('uint16')
        elif self.nbits == 8:
            out = block.flatten(order='C').astype('uint8')
        elif self.nbits in [1, 2, 4]:
            raw = block.flatten(order='C').astype('uint8')
            values_per_byte = 8 // self.nbits

            trimmed_len = (len(raw) // values_per_byte) * values_per_byte
            raw = raw[:trimmed_len]
            
            reshaped_raw = raw.reshape(-1, values_per_byte)
            out = np.zeros(len(reshaped_raw), dtype=np.uint8)

            for i in range(values_per_byte):
                out |= np.left_shift(reshaped_raw[:, i], i * self.nbits)
        else:
            sys.exit(f"{self.nbits} bit filterbank is not supported.")

        self.write_file.write(out.tobytes())


def merge_filterbanks(filterbanks, output_file, gulp_size=2**11):
    fb_list = [FilterbankReader(fb, stats_samples=0) for fb in filterbanks]
    writer = FilterbankWriter(fb_list[0], output_file)

    for fb in fb_list:
        n_blocks, remainder = divmod(fb.n_samples, gulp_size)

        for _ in range(n_blocks):
            read_block = fb.read_block(gulp_size)
            writer.write_block(read_block)

        read_block = fb.read_block(remainder)
        writer.write_block(read_block)

    del fb_list, writer


def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + str(output))

def read_datfile(path, nbits, count=-1):
    if nbits == 64:
        raw = np.fromfile(path, dtype=np.float64,count=count)
        return raw.astype(np.float64)
    if nbits == 32:
        raw = np.fromfile(path, dtype=np.float32,count=count)
        return raw.astype(np.float64)
    if nbits == 8:
        raw = np.fromfile(path, dtype=np.uint8, count=count)
        return raw.astype(np.uint8)