import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


def profile(t, phase, sigma, Amp):
    g = norm(phase, sigma).pdf(t)
    return Amp * g

def profile_2(t, phase1, sigma1, Amp1, phase2, sigma2, Amp2):
    g1 = norm(phase1, sigma1).pdf(t)
    g2 = norm(phase2, sigma2).pdf(t)
    return Amp1 * g1 + Amp2 * g2


def scale_freq_phase(freq_phase, intensity_profile):
    
    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile)

    def fit_phase(phase, phase_off, Amp1, Amp2, base):
        g = profile_2((phase-phase_off)%1, *profile_pars[:2], Amp1, *profile_pars[3:5], Amp2)
        return g + base

    prof_nbins = len(intensity_profile)
    nchans = len(freq_phase.T)
    phase_arr = np.linspace(0, 1, prof_nbins)
    freq_phase = np.roll(freq_phase, prof_nbins//2-phase_corr, axis=1)

    profile_s = []
    for i in np.arange(nchans):
        time_phase_arr = freq_phase[i]  

        out = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 1, 0], bounds =[[-profile_pars[1], 0, 0, -np.inf], [profile_pars[1], np.inf, np.inf, np.inf]])
        
        updated_pars = profile_pars.copy()
        updated_pars[2] = out[0][1]
        updated_pars[5] = out[0][2]
        spectra = profile_2(phase_arr-out[0][0], *updated_pars)
        profile_s.append(spectra)

    return np.array(profile_s)


def get_IP_interp(intensity_profile):

    true_prof = intensity_profile
    true_prof /= np.max(true_prof)

    max_ind = np.argmax(true_prof)

    true_prof = np.roll(true_prof, len(true_prof)//2-max_ind)
    phase = np.linspace(0, 1, len(true_prof))

    out = curve_fit(profile_2, phase, true_prof, p0=[0.49, 0.05, 0.5, 0.51, 0.05, 0.5], 
                    bounds=[[0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1]])

    noise = true_prof-profile_2(phase, *out[0])
    SNR = np.sqrt(np.sum((true_prof-np.mean(noise))**2)/np.std(noise)**2)

    return out[0], max_ind, SNR


def fit_time_phase(time_phase, intensity_profile, obs_len):

    profile_pars, phase_corr, SNR = get_IP_interp(intensity_profile)

    def fit_phase(phase, phase_off, Amp, base):
        g = profile_2((phase-phase_off)%1, *profile_pars[:2], profile_pars[2]*Amp, *profile_pars[3:5], profile_pars[5]*Amp)
        return g + base

    prof_nbins = len(intensity_profile)
    time_nbins = len(time_phase)
    phase_arr = np.linspace(0, 1, prof_nbins)
    time_phase /= np.max(time_phase)
    time_phase = np.roll(time_phase, prof_nbins//2-phase_corr, axis=1)

    time_amp, phase, err = [], [], []
    for i in np.arange(time_nbins):
        time_phase_arr = time_phase[i]

        out = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 0])
        phase.append((out[0][0]+0.5)%1-0.5)
        time_amp.append(out[0][1])

        fit_err =np.sqrt(np.diag(out[1]))[0]
        subint_err = (np.std(time_phase_arr-fit_phase(phase_arr, *out[0]))/(1-out[0][2]))**2
        err.append(np.sqrt(fit_err**2 + subint_err**2))


    def fit_f(t, theta0, f0, f1, f2, f3):
        return theta0 + f0*t + (1/2)*f1*t**2 + (1/6)*f2*t**3 + (1/24)*f3*t**4

    dt = obs_len / time_nbins
    time = (np.arange(time_nbins) + 0.5) * dt - obs_len/2
    out = curve_fit(fit_f, time, phase, sigma=np.array(err), p0=[1e-3,1e-6,1e-8,1e-10,1e-12])
    phase_offset = -out[0][0] #- (phase_arr[phase_corr] - phase_arr[prof_nbins//2])
    phase_shift = (phase_arr[phase_corr] - phase_arr[prof_nbins//2])

    freq_deriv = {}
    for i, fx in enumerate(out[0][1:]):
        freq_deriv[f'F{i}'] = -fx

    time_amp_smooth = savgol_filter(time_amp, window_length=11, polyorder=3)

    return freq_deriv, phase_offset, phase_shift, SNR, time+obs_len/2, time_amp_smooth


def fit_phase_offset(intensity_profile_OPT, intensity_profile_INIT):

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)

    prof_nbins = len(intensity_profile_OPT)
    phase = np.linspace(0, 1, prof_nbins)
    intensity_profile_OPT -= np.min(intensity_profile_OPT)
    intensity_profile_OPT /= np.max(intensity_profile_OPT)
    intensity_profile_OPT = np.roll(intensity_profile_OPT, prof_nbins//2-phase_corr)

    def func_r(x, A1, A2, x1, x2, d):
        g1 = profile_2((x-x1) % 1, *profile_pars)
        g2 = profile_2((x-x2) % 1, *profile_pars)
        return A1*g1 + A2*g2 + d
    
    p0 = [0.5, -0.5, 0, phase[np.argmin(intensity_profile_OPT)]-phase[np.argmax(intensity_profile_OPT)], 0.5]
    out = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0)

    phase_offset = out[0][3]-out[0][2]
    SNR_scale = np.abs(out[0][0]/out[0][1])
    
    return phase_offset, SNR_scale
    

