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
    from scipy.stats import norm
    from scipy.optimize import curve_fit

    profile_pars, phase_corr, SNR = get_IP_interp(intensity_profile)

    def profile(t, phase, sigma, Amp):
        g = norm(phase, sigma).pdf(t)
        return Amp * g/g.max()

    def fit_phase(phase, phase_off, Amp, base):
        g = profile((phase-phase_off)%1, *profile_pars[:2], Amp)
        return g + base

    prof_nbins = len(intensity_profile)
    nchans = len(freq_phase.T)
    phase_arr = np.linspace(0, 1, prof_nbins)
    freq_phase = np.roll(freq_phase, prof_nbins//2-phase_corr, axis=1)

    profile_s = []
    for i in np.arange(nchans):
        time_phase_arr = freq_phase[i]

        out = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 0], bounds =[[-profile_pars[1], 0, -np.inf], 
                                                                                        [profile_pars[1], np.inf, np.inf]])
        
        spectra = out[0][1]*profile(phase_arr-out[0][0], *profile_pars)
        profile_s.append(spectra)

    return np.array(profile_s)

def get_IP_interp(intensity_profile):
    import numpy as np
    from scipy.stats import norm
    from scipy.optimize import curve_fit

    def profile(t, phase, sigma, Amp):
        g = norm(phase, sigma).pdf(t)
        return Amp * g/g.max()

    true_prof = intensity_profile
    true_prof /= np.max(true_prof)

    max_ind = np.argmax(true_prof)

    true_prof = np.roll(true_prof, len(true_prof)//2-max_ind)
    phase = np.linspace(0, 1, len(true_prof))

    out = curve_fit(profile, phase, true_prof)

    noise = true_prof-profile(phase, *out[0])
    SNR = np.sqrt(np.sum((true_prof-np.mean(noise))**2)/np.std(noise)**2)

    return out[0], max_ind, SNR


def fit_time_phase(time_phase, intensity_profile, obs_len):
    import numpy as np
    from scipy.stats import norm
    from scipy.optimize import curve_fit

    profile_pars, phase_corr, SNR = get_IP_interp(intensity_profile)

    def profile(t, phase, sigma, Amp):
        g = norm(phase, sigma).pdf(t)
        return Amp * g/g.max()

    def fit_phase(phase, phase_off, Amp, base):
        g = profile((phase-phase_off)%1, *profile_pars[:2], Amp)
        return g + base

    prof_nbins = len(intensity_profile)
    time_nbins = len(time_phase)
    phase_arr = np.linspace(0, 1, prof_nbins)
    time_phase = np.roll(time_phase, prof_nbins//2-phase_corr, axis=1)

    phase, err = [], []
    for i in np.arange(time_nbins):
        time_phase_arr = time_phase[i]

        out = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 0])
        phase.append((out[0][0]+0.5)%1-0.5)

        fit_err =np.sqrt(np.diag(out[1]))[0]
        subint_err = (np.std(time_phase_arr-fit_phase(phase_arr, *out[0]))/(1-out[0][2]))**2
        err.append(np.sqrt(fit_err**2 + subint_err**2))


    def fit_f(t, theta0, f0, f1, f2, f3):
        return theta0 + f0*t + (1/2)*f1*t**2 + (1/6)*f2*t**3 + (1/24)*f3*t**4

    time = np.linspace(-obs_len/2, obs_len/2, time_nbins)
    out = curve_fit(fit_f, time, phase, sigma=np.array(err), p0=[1e-3,1e-6,1e-8,1e-10,1e-12])
    phase_offset = -out[0][0] - (phase_arr[phase_corr] - phase_arr[prof_nbins//2])

    freq_deriv = {}
    for i, fx in enumerate(out[0][1:]):
        freq_deriv[f'F{i}'] = -fx

    return freq_deriv, phase_offset, SNR


def fit_phase_offset(intensity_profile_OPT, intensity_profile_INIT):
    import numpy as np
    from scipy.stats import norm
    from scipy.optimize import curve_fit

    profile_pars, phase_corr, SNR = get_IP_interp(intensity_profile_INIT)

    def profile(t, phase, sigma, Amp):
        g = norm(phase, sigma).pdf(t)
        return Amp * g/g.max()

    prof_nbins = len(intensity_profile_OPT)
    phase = np.linspace(0, 1, prof_nbins)
    intensity_profile_OPT -= np.min(intensity_profile_OPT)
    intensity_profile_OPT /= np.max(intensity_profile_OPT)
    intensity_profile_OPT = np.roll(intensity_profile_OPT, prof_nbins//2-phase_corr, axis=1)

    def func_r(x, A1, A2, x1, x2, d):
        g1 = profile((x-x1) % 1, *profile_pars)
        g2 = profile((x-x2) % 1, *profile_pars)
        return A1*g1/g1.max() + A2*g2/g2.max() + d
    
    p0 = [0.5, -0.5, 0, phase[np.argmin(intensity_profile_OPT)]-phase[np.argmax(intensity_profile_OPT)], 0.5]
    out = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0)

    phase_offset = out[0][3]-out[0][2]
    SNR_scale = np.abs(out[0][0]/out[0][1])
    
    return phase_offset, SNR_scale
    

