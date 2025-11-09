# modified verion of:

# feature extractor Copyright (C) 2025  Devika Bhatnagar

# This program comes with ABSOLUTELY NO WARRANTY;
# This is free software, and you are welcome to redistribute it under certain conditions;
# https://github.com/erc-compact/new_feature_extractor


import re, os
import subprocess
import numpy as np


class ARProcessor:
    def __init__(self, file, mode='extract', work_dir='cwd', search=False):
        self.ar_file = file
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        if mode == 'extract':
            self.extract(search)
        else:
            self.load(file)

    def extract(self, search):
        self.run_dmffdot(search)
        self.fits_file = self.ar_file.replace('.ar', '.px')

    def load(self, fits_path):
        from astropy.io import fits

        self.hdul = fits.open(fits_path)

    def get_time_phase(self):
        time_phase = self.hdul[11].data[0]
        x = np.array(time_phase[0])
        y = np.array(time_phase[1])
        z = np.array(time_phase[2])
        TP = z.reshape((y.shape[0],x.shape[0]))
        return TP[:, :self.TP.shape[1]//2] 

    def get_freq_phase(self):
        freq_phase = self.hdul[8].data[0]
        x = np.array(freq_phase[0])
        y = np.array(freq_phase[1])
        z = np.array(freq_phase[2])
        FP = z.reshape((y.shape[0],x.shape[0])) 
        return FP[:, :self.FP.shape[1]//2] 

    def get_dm_curve(self):
        self.dm_curve = self.hdul[21].data[0][1] 
    
    def get_intensity_prof(self):
        intensity_prof = self.hdul[7].data[0][1] 
        return intensity_prof[:len(self.intensity_prof)//2]

    def get_ffdot(self):
        ffdot = self.hdul[16].data[0]
        x = np.array(ffdot[0])
        y = np.array(ffdot[1])
        z = np.array(ffdot[2])
        self.ffdot = z.reshape((x.shape[0],y.shape[0]))

    @staticmethod
    def extract_number(text):
        match = re.search(r"[-+]?\d*\.\d+|\d+", text.split('=')[1])  
        return match.group(0) if match else None

    def extract_and_convert(self, keyword):
        value = next((self.extract_number(row[4]) for row in self.hdul[15].data if keyword in row[4]), None)
        return float(value) if value is not None else None

    def get_period(self):
        return self.extract_and_convert("P0")

    def get_pepoch(self):
        return self.extract_and_convert("Pepoch")

    def get_DM(self):
        return self.extract_and_convert("DM")

    def get_SNR(self):
        return self.extract_and_convert("S/N")
    

    def run_dmffdot(self, search):
        if search:
            command = ["dmffdot", "-f", self.ar_file, "--saveimage"]
        else:
            command = ["dmffdot", "-f", self.ar_file, "--saveimage", "--nosearch"]

        print("Running command:", " ".join(command))

        subprocess.run(command, cwd=self.work_dir) 

        print("Run complete")
        

