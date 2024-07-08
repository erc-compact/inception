import struct
import numpy as np 


class FilterbankIO:
    """ Edited version of Mike Keith's filtool code"""
    _inttypes = ['machine_id', 'telescope_id', 'data_type', 'nchans',\
                    'nbits', 'nifs', 'scan_number', 'barycentric','pulsarcentric',\
                    'nbeams', 'ibeam', 'nsamples']
    _strtypes = ['source_name','rawdatafile']
    _dbltypes = ['tstart','tsamp','fch1','foff','refdm','az_start','za_start','src_raj','src_dej']
    _chrtypes = ['signed']

    def __init__(self, filterbank):
        self.header={}
        self.read_file = None
        self.write_file = None
        self.read_header_pos = 0
        self.read_data_pos = 0
        self.write_header_pos = 0
        self.write_data_pos = 0

        self.load(filterbank)
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
            while key != "HEADER_END":
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
    
    def get_n_samples(self):
        self.read_file.seek(0,2)
        nchans = self.header['nchans']
        nbits = self.header['nbits']
        n_samples = int((self.read_file.tell() - self.read_data_pos)/nchans * (8/nbits))
        self.read_file.seek(self.read_data_pos)
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

    def read_block(self, nsamples, nbits=None):
        if nbits == None:
            nbits = self.header['nbits']
        nchans = self.header['nchans']
        dtype=np.float64

        if nbits == 64:
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file,dtype=np.float64,count=nbytes)
            return raw.reshape(-1,nchans).astype(dtype)
        if nbits == 8:
            nbytes = nsamples * nchans
            raw = np.fromfile(self.read_file, dtype=np.uint8, count=nbytes)
            return raw.reshape(-1, nchans).astype(dtype)
        
    def write_block(self, block, nbits=None): 
        if nbits == None:
            nbits = self.header['nbits']
            
        if nbits == 64:
            raw = block.flatten(order='C').astype('float64')
            self.write_file.write(raw.tobytes())
        elif nbits == 8:
            max_val = 2**nbits-1
            block[block>max_val] = max_val
            block[block<0] = 0
            raw = block.flatten(order='C').astype('uint8')
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
    
    def write_channel(self, channel_number, input_path, nbits=None):
        if nbits == None:
            nbits = self.header['nbits']
        channel_data = read_datfile(input_path, nbits)
        number_of_channels = self.header["nchans"]
        for byte_i in range(len(channel_data)):
            offset = byte_i*number_of_channels + (channel_number%number_of_channels)
            self.read_file.seek(self.read_data_pos + offset)
            self.read_file.write(channel_data[byte_i].tobytes())
        self.read_file.seek(self.read_data_pos)


def read_datfile(path, nbits, count=-1):
    if nbits == 64:
        raw = np.fromfile(path, dtype=np.float64,count=count)
        return raw.astype(np.float64)
    if nbits == 8:
        raw = np.fromfile(path, dtype=np.uint8, count=count)
        return raw.astype(np.uint8)