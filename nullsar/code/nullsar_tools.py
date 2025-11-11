import sys, os
import json
import subprocess

def parse_JSON(json_file):
    try:
        with open(json_file, 'r') as file:
            pars = json.load(file)
    except FileNotFoundError:
        sys.exit(f'Unable to find {json_file}.')
    except json.JSONDecodeError:
        sys.exit(f'Unable to parse {json_file} using JSON.')
    else:
        return pars
    

def execute(cmd):
    os.system(cmd)

def print_exe(output):
    execute("echo " + str(output))
    

def rsync(source, destination, shell=True):
    try:
        subprocess.run(f'rsync -Pav {source} {destination}', shell=shell, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        subprocess.run(f'cp -av {source} {destination}', shell=shell)


def parse_par_file(parfile):
    data = {}
    with open(parfile, "r") as f:
        for line in f:
            columns = line.split()
            data[columns[0]] = columns[1]
    return data


def parse_cand_file(candifle):
    with open(candifle, 'r') as f:
        data = [line.strip().split() for line in f if line.strip()]
    return dict(zip(data[0], data[1]))


def scale_freq_phase(freq_phase, intensity_profile):
    import numpy as np

    mask = (intensity_profile < np.max(intensity_profile)*0.05)

    freq_phase.T[mask] = 0
    freq_phase /= np.max(freq_phase)
    freq_phase[freq_phase < 0] = 0

    return freq_phase

def get_IP_interp(freq_phase):
    import numpy as np
    from scipy.interpolate import interp1d

    true_prof = np.sum(freq_phase, axis=0)
    true_prof -= np.min(true_prof)
    true_prof /= np.max(true_prof)

    phase = np.linspace(0, 1, len(true_prof))
    return phase, interp1d(phase, true_prof, kind='quadratic')


def fit_time_phase(time_phase, freq_phase, obs_len):
    import numpy as np
    from scipy.optimize import curve_fit

    phase_arr, profile = get_IP_interp(freq_phase)

    def func(t, phase_off, Amp, base):
        g = profile((t-phase_off)%1)
        return Amp * g/g.max() + base

    phase, err = [], []
    for i in np.arange(len(time_phase)):
        time_phase_arr = time_phase[i]
        time_phase_arr -= np.min(time_phase_arr)
        time_phase_arr /= np.max(time_phase_arr)

        out = curve_fit(func, phase_arr, time_phase_arr, p0=[0, np.max(time_phase_arr), np.min(time_phase_arr)])
        phase.append(out[0][0])
        err.append(np.sqrt(np.diag(out[1]))[0])

    def fit_f(t, theta0, f0, f1, f2, f3):
        return theta0 + f0*t + (1/2)*f1*t**2 + (1/6)*f2*t**3 + (1/24)*f3*t**4

    time = np.linspace(-obs_len/2, obs_len/2, len(time_phase))
    out = curve_fit(fit_f, time, phase, sigma=np.array(err), p0=[1e-3,1e-6,1e-8,1e-10,1e-12])
    phase_offset = -out[0][0]

    freq_deriv = {}
    for i, fx in enumerate(out[0][1:]):
        freq_deriv[f'F{i}'] = -fx

    return freq_deriv, phase_offset


def fit_phase_offset(intensity_profile, freq_phase):
    import numpy as np
    from scipy.optimize import curve_fit

    phase, profile = get_IP_interp(freq_phase)

    intensity_profile -= np.min(intensity_profile)
    intensity_profile /= np.max(intensity_profile)

    def func_r(x, A1, A2, x2, d):
        g1 = profile(x % 1)
        g2 = profile((x-x2) % 1)
        return A1* g1/g1.max() +A2* g2/g2.max()+ d
    
    p0 = [0.5, -0.5, phase[np.argmin(intensity_profile)]-phase[np.argmax(intensity_profile)], 0.5]
    out = curve_fit(func_r, phase, intensity_profile, p0=p0)

    phase_offset = out[0][2]
    SNR_scale = np.abs(out[0][0]/out[0][1])
    
    return phase_offset, SNR_scale
    

