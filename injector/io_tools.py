import os
import sys
import struct
import numpy as np 
# update and add io_tools to new folder in inception + extract_scripts

class FilterbankIO: # edited version from Michael Keith's filtools
    _inttypes = ['machine_id', 'telescope_id', 'data_type', 'nchans',\
                    'nbits', 'nifs', 'scan_number', 'barycentric','pulsarcentric',\
                    'nbeams', 'ibeam', 'nsamples']
    _strtypes = ['source_name','rawdatafile']
    _dbltypes = ['tstart','tsamp','fch1','foff','refdm','az_start','za_start','src_raj','src_dej']
    _chrtypes = ['signed']

    def __init__(self, filterbank):
        self.header={}
        self.path = filterbank
        self.read_file = None
        self.write_file = None
        self.read_header_pos = 0
        self.read_data_pos = 0
        self.write_header_pos = 0
        self.write_data_pos = 0

        self.load(filterbank)
        self.nchans = self.header['nchans']
        self.fb_mean, self.fb_std = self.get_FB_stats(2**14)
        self.n_samples = self.get_n_samples()

    def load(self, filename):
        def read_string(file):
            nchar = np.fromfile(file,dtype=np.int32,count=1)[0]
            if nchar == 0:
                return ""
            if nchar < 1 or nchar > 80:
                raise Exception(f"Cannot parse filterbank header (Nchar was {nchar} when reading string).")
            byte_data = file.read(nchar)
            string_data = byte_data.decode("UTF-8")
            return string_data
        
        self.read_file = open(filename, 'r+b')
        self.read_header_pos=self.read_file.tell()

        key = read_string(self.read_file)
        bytes_read = len(key)+4
        if key=="HEADER_START":
            key = read_string(self.read_file)
            bytes_read += len(key)+4
            while key != "HEADER_END": # optimise
                if key in FilterbankIO._inttypes:
                    self.header[key] = np.fromfile(self.read_file,dtype=np.int32,count=1)[0]
                    bytes_read += 4
                elif key in FilterbankIO._strtypes:
                    self.header[key] = read_string(self.read_file)
                    bytes_read += len(self.header[key])+4
                elif key in FilterbankIO._dbltypes:
                    self.header[key] = np.fromfile(self.read_file,dtype=np.float64,count=1)[0]
                    bytes_read += 8
                elif key in FilterbankIO._chrtypes:
                    self.header[key] = np.fromfile(self.read_file,dtype=np.int8,count=1)[0]
                    bytes_read += 1
                else:
                    raise Exception(f"Cannot parse filterbank header, key '{key}' not understood")
                key = read_string(self.read_file)
                bytes_read += len(key)+4
        else:
            raise Exception("Cannot parse filterbank header, HEADER_START was not found")
        self.read_data_pos=self.read_file.tell()
        
        return bytes_read

    def get_FB_stats(self, nsamples):
        cuurent_loc = self.read_file.tell()
        data = self.read_block(nsamples)
        self.read_file.seek(cuurent_loc, 0)
        return np.median(np.mean(data, axis=0)), np.median(np.std(data, axis=0))
    
    def get_n_samples(self, which_file='read'):
        file = self.read_file if which_file == 'read' else self.write_file
        file_dat = self.read_data_pos if which_file == 'read' else self.write_data_pos
        file.seek(0,2)
        nchans = self.header['nchans']
        nbits = self.header['nbits']
        n_samples = int((file.tell() - file_dat)/nchans * (8/nbits))
        self.read_file.seek(file_dat)
        return n_samples
    
    def new_filterbank(self, filename, nbits=None):
        if nbits == None:
            nbits = self.header['nbits']
        else:
            self.header['nbits'] = nbits

        def write_string(file, string):
            l=len(string)
            byte_data= struct.pack(f"i{l}s", l, string.encode("UTF-8"))
            file.write(byte_data)
        
        self.write_file = open(filename, "wb")
        self.write_header_pos=self.write_file.tell()

        write_string(self.write_file,"HEADER_START")
        for key in self.header:
            write_string(self.write_file,key)
            if key in FilterbankIO._inttypes:
                self.write_file.write(struct.pack("i",self.header[key]))
            elif key in FilterbankIO._strtypes:
                write_string(self.write_file,self.header[key])
            elif key in FilterbankIO._dbltypes:
                self.write_file.write(struct.pack("d",self.header[key]))
            elif key in FilterbankIO._chrtypes:
                self.write_file.write(struct.pack("b",self.header[key]))
            else:
                raise Exception(f"Cannot understand filterbank header for writing, key '{key}' not understood")
        write_string(self.write_file,"HEADER_END")
        self.write_data_pos=self.write_file.tell()

    def get_reader(self):
        nchans = self.nchans
        dtype=np.float64
        
        def read_block_32bit(nsamples,dtype=dtype):
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file,dtype=np.float32,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        
        def read_block_16bit(nsamples,dtype=dtype):
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file,dtype=np.uint16,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        
        def read_block_8bit(nsamples,dtype=dtype):
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file,dtype=np.uint8,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        
        def read_block_4bit(nsamples,dtype=dtype):
            nbytes = int((nsamples * nchans)/2)
            raw = np.fromfile(self.read_file,dtype=np.uint8,count=nbytes)    
            out = np.zeros(len(raw) * 2, dtype=dtype)
            out[::2] += np.bitwise_and(raw,0x0f)
            out[1::2] += np.right_shift(raw,4)
            return out.reshape(-1,nchans).astype(dtype)       
        
        def read_block_2bit(nsamples,dtype=dtype):
            nbytes = int((nsamples * nchans)/4)
            raw = np.fromfile(self.read_file,dtype=np.uint8,count=nbytes)
            out = np.zeros(len(raw) * 4, dtype=dtype)
            out[::4] += np.bitwise_and(raw,0x03)
            out[1::4] += np.right_shift(np.bitwise_and(raw,0x0f),2)
            out[2::4] += np.right_shift(np.bitwise_and(raw,0x3f),4)
            out[3::4] += np.right_shift(np.bitwise_and(raw,0xff),6)
            return out.reshape(-1,nchans).astype(dtype)
        
        def read_block_1bit(nsamples,dtype=dtype):
            nbytes = int((nsamples * nchans)/8)
            raw = np.fromfile(self.read_file,dtype=np.uint8,count=nbytes)
            
            out = np.zeros(len(raw) * 8, dtype=dtype)
            out[::8] += np.bitwise_and(raw,0x01)
            out[1::8] += np.right_shift(np.bitwise_and(raw,0x03),1)
            out[2::8] += np.right_shift(np.bitwise_and(raw,0x07),2)
            out[3::8] += np.right_shift(np.bitwise_and(raw,0x0f),3)
            out[4::8] += np.right_shift(np.bitwise_and(raw,0x1f),4)
            out[5::8] += np.right_shift(np.bitwise_and(raw,0x3f),5)
            out[6::8] += np.right_shift(np.bitwise_and(raw,0x7f),6)
            out[7::8] += np.right_shift(np.bitwise_and(raw,0xff),7)
            return out.reshape(-1,nchans).astype(dtype)
        
        
        if 'signed' in self.header and self.header['signed']!= 1:
            sys.exit(f'Filterbank data must be unsigned')
        elif self.header['nbits'] == 16:
            return read_block_16bit
        elif self.header['nbits'] == 8:
            return read_block_8bit
        elif self.header['nbits'] == 4:
            return read_block_4bit
        elif self.header['nbits'] == 2:
            return read_block_2bit
        elif self.header['nbits'] == 1:
            return read_block_1bit
        else:
            sys.exit(f"{self.header['nbits']} bit filterbank is not supported.")

    def get_writer(self):
        
        def write_block_32bit(block):
            raw = block.flatten(order='C').astype('float32')
            self.write_file.write(raw.tobytes())
            
        def write_block_16bit(block):
            block[block>65535] = 65535
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint16')

            self.write_file.write(raw.tobytes())
            
        def write_block_8bit(block):
            block[block>255] = 255
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint8')
            
            self.write_file.write(raw.tobytes())

        def write_block_4bit(block):
            block[block > 15] = 3
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint8')
            
            out = np.zeros(int(len(raw)/2),dtype=np.uint8)
            out += raw[::2]
            out += np.left_shift(raw[1::2],4)
            self.write_file.write(out.tobytes())
            
        def write_block_2bit(block):
            block[block > 3] = 3
            block[block < 0] = 0
            
            raw = block.flatten(order='C').astype('uint8')
            out = np.zeros(int(len(raw)/4),dtype=np.uint8)
            out += raw[::4]
            out += np.left_shift(raw[1::4],2)
            out += np.left_shift(raw[2::4],4)
            out += np.left_shift(raw[3::4],6)
            self.write_file.write(out.tobytes())
            
        def write_block_1bit(block):
            block[block > 1] = 1
            block[block < 0] = 0
            raw = block.flatten(order='C').astype('uint8')
            out = np.zeros(int(len(raw)/8),dtype=np.uint8)
            out += raw[::8]
            out += np.left_shift(raw[1::8],1)
            out += np.left_shift(raw[2::8],2)
            out += np.left_shift(raw[3::8],3)
            out += np.left_shift(raw[4::8],4)
            out += np.left_shift(raw[5::8],5)
            out += np.left_shift(raw[6::8],6)
            out += np.left_shift(raw[7::8],7)
            self.write_file.write(out.tobytes())

            
        if 'signed' in self.header and self.header['signed']!= 1:
            sys.exit(f'Filterbank data must be unsigned')     
        elif self.header['nbits'] == 16:
            return write_block_16bit
        elif self.header['nbits'] == 8:
            return write_block_8bit
        elif self.header['nbits'] == 4:
            return write_block_4bit
        elif self.header['nbits'] == 2:
            return write_block_2bit
        elif self.header['nbits'] == 1:
            return write_block_1bit
        else:
            sys.exit(f"{self.header['nbits']} bit filterbank is not supported.")

    def read_block(self, nsamples):   
        return self.get_reader()(nsamples, dtype=np.float64)
    
    def write_block(self, block):   
        return self.get_writer()(block)

    def read_block_old(self, nsamples, nbits=None):
        if nbits == None:
            nbits = self.header['nbits']
        nchans = self.header['nchans']
        dtype=np.float64

        if nbits == 64:
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file,dtype=np.float64,count=nbytes)
            return raw.reshape(-1, nchans).astype(dtype)
        if nbits == 8:
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file, dtype=np.uint8, count=nbytes)
            return raw.reshape(-1, nchans).astype(dtype)
        
    def write_block_old(self, block, nbits=None): 
        if nbits == None:
            nbits = self.header['nbits']  
        if nbits == 64:
            raw = block.flatten(order='C').astype(np.float64)
            self.write_file.write(raw.tobytes())
        elif nbits == 8:
            max_val = 2**nbits-1
            block[block>max_val] = max_val
            block[block<0] = 0
            raw = block.flatten(order='C').astype(np.uint8)
            self.write_file.write(raw.tobytes())

    def split_block(self, output_path, block, mode, nbits=None):
        if nbits == None:
            nbits = self.header['nbits']
        nchans = self.header['nchans']
        for chan in range(nchans):
            with open(output_path+f'_{chan}.dat', mode) as file:
                if nbits == 64:
                    raw = block.T[chan].astype('float64')
                elif nbits == 8:
                    raw = block.T[chan].astype('uint8')
                file.write(raw.tobytes())
    
    def split_fb_channels(self, output_path, ext=''):
        nchans = self.nchans

        block_size = 2**11
        n_blocks, remainder = divmod(self.n_samples, block_size)
        open_files = [open(output_path+f'/{ext}_{chan}.dat', 'wb') for chan in range(nchans)]

        def split_block(block):
            for chan in range(nchans):
                raw = block.T[chan].astype('float64')
                open_files[chan].write(raw.tobytes())

        for _ in range(n_blocks):
            block = self.read_block(block_size)
            split_block(block)
        
        block_remain = self.read_block(remainder)
        split_block(block_remain)

def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + str(output))

def read_datfile(path, nbits, count=-1):
    if nbits == 64:
        raw = np.fromfile(path, dtype=np.float64,count=count)
        return raw.astype(np.float64)
    if nbits == 8:
        raw = np.fromfile(path, dtype=np.uint8, count=count)
        return raw.astype(np.uint8)

def str2func(value, par, id, func):
    try:
        converted_value = func(value)
    except ValueError:
        sys.exit(f'Error: Invalid {par} for pulsar {id}')
    return converted_value


