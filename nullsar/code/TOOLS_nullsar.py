import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.interpolate import PchipInterpolator



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
    amps = []
    for i in np.arange(nchans):
        time_phase_arr = freq_phase[i]  
        try:
            out = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 1, 0], bounds =[[-profile_pars[1], 0, 0, -np.inf], [profile_pars[1], np.inf, np.inf, np.inf]])
        except:
            amps.append((0, 0))
            spectra = profile_2(phase_arr, *profile_pars)
            profile_s.append(spectra)
        else:
            updated_pars = profile_pars.copy()
            updated_pars[2] = out[0][1]
            updated_pars[5] = out[0][2]
            amps.append((out[0][1], out[0][2]))
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
                    bounds=[[0, 0.01, 0, 0, 0.01, 0], [1, 1, 1, 1, 1, 1]])

    noise = true_prof-profile_2(phase, *out[0])
    SNR = np.sqrt(np.sum((true_prof-np.mean(noise))**2)/np.std(noise)**2)

    return out[0], max_ind, SNR


def fit_f(t, theta0, f0, f1):
    return theta0 + f0*t + (1/2)*f1*t**2


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

    dt = obs_len / time_nbins
    time = (np.arange(time_nbins) + 0.5) * dt - obs_len/2
    out = curve_fit(fit_f, time, phase, sigma=np.array(err), p0=[1e-3,1e-6,1e-10])
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
    
    p0 = [0.5, -0.5, 0.5-phase[np.argmax(intensity_profile_OPT)], phase[np.argmin(intensity_profile_OPT)]-0.5, 0.5]
    out = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0)

    phase_offset = out[0][3]-out[0][2]
    SNR_scale = np.abs(out[0][0]/out[0][1])
    
    return phase_offset, SNR_scale, out


def fit_subint_phase_offset(time_phase_OPT, intensity_profile_INIT, phase_off):

    def fit_subint(intensity_profile_OPT, func_r, phase_off):
        prof_nbins = len(intensity_profile_OPT)
        phase = np.linspace(0, 1, prof_nbins)

        if np.std(intensity_profile_OPT) == 0:
            return 1

        intensity_profile_OPT -= np.min(intensity_profile_OPT)
        intensity_profile_OPT /= np.max(intensity_profile_OPT)
        intensity_profile_OPT = np.roll(intensity_profile_OPT, prof_nbins//2-np.argmax(intensity_profile_OPT))
        
        p0 = [0.5, -0.5, 0, phase_off, np.median(intensity_profile_OPT)]
        bounds = [[0, -1, -0.4, -0.5, 0],[1, 0, 0.4, 0.5, 1]]
        try:
            out = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0, bounds=bounds)
        except:
            return 1
        else:
            phase_offset = out[0][3]-out[0][2]
            SNR_scale = np.abs(out[0][0]/out[0][1])

            SNR_A = abs(out[0][0]/np.std(intensity_profile_OPT-func_r(phase, *out[0])))
            SNR_B = abs(out[0][1]/np.std(intensity_profile_OPT-func_r(phase, *out[0])))

            if (SNR_A  > 5) & (SNR_B > 5):
                return SNR_scale
            else:
                return 1

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)
    
    def func_r(x, A1, A2, x1, x2, d):
        g1 = profile_2((x-x1) % 1, *profile_pars)
        g2 = profile_2((x-x2) % 1, *profile_pars)
        return A1*g1 + A2*g2 + d

    snr_arr = []
    for i, subint_i in enumerate(time_phase_OPT):
        SNR_scale = fit_subint(subint_i, func_r, phase_off)
        snr_arr.append(SNR_scale)
    
    time_amp_smooth = savgol_filter(snr_arr, window_length=11, polyorder=3)
    return time_amp_smooth


def fit_chan_phase_offset(freq_phase_OPT, intensity_profile_INIT, freq_arr, phase_off, PSR_P0):

    def fit_chan(intensity_profile_OPT, func_r, phase_off):
        prof_nbins = len(intensity_profile_OPT)
        phase = np.linspace(0, 1, prof_nbins)

        if np.std(intensity_profile_OPT) == 0:
            return phase_off, 1000

        intensity_profile_OPT -= np.min(intensity_profile_OPT)
        intensity_profile_OPT /= np.max(intensity_profile_OPT)
        intensity_profile_OPT = np.roll(intensity_profile_OPT, prof_nbins//2-np.argmax(intensity_profile_OPT))
        
        p0 = [0.5, -0.5, 0, phase_off, np.median(intensity_profile_OPT)]
        bounds = [[0, -1, -0.4, -0.5, 0],[1, 0, 0.4, 0.5, 1]]
        try:
            out = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0, bounds=bounds)
        except:
            return phase_off, 1000
        else:
            phase_offset = out[0][3]-out[0][2]
            SNR_scale = np.abs(out[0][0]/out[0][1])

            SNR_A = abs(out[0][0]/np.std(intensity_profile_OPT-func_r(phase, *out[0])))
            SNR_B = abs(out[0][1]/np.std(intensity_profile_OPT-func_r(phase, *out[0])))

            if (SNR_A  > 5) & (SNR_B > 5):
                return phase_offset, (2/(SNR_A+SNR_B))**2
            else:
                return phase_off, 1000

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)
    
    def func_r(x, A1, A2, x1, x2, d):
        g1 = profile_2((x-x1) % 1, *profile_pars)
        g2 = profile_2((x-x2) % 1, *profile_pars)
        return A1*g1 + A2*g2 + d

    def delay(freq_arr, DM, off):
        DM_const = 4148.80642389
        return -DM * DM_const / freq_arr**2 + off

    phase_arr = []
    phase_sigma = []
    for i, chan_i in enumerate(freq_phase_OPT):
        phase_offset, phase_err = fit_chan(chan_i, func_r, phase_off)
        phase_arr.append(phase_offset)
        phase_sigma.append(phase_err)
    
    delays = np.array(phase_arr) * PSR_P0
    dm_out = curve_fit(delay, freq_arr, delays, sigma=phase_sigma)

    sig_level = dm_out[0][0]/np.sqrt(np.diag(dm_out[1]))[0]
    if sig_level > 3:
        DM_offset = dm_out[0][0]
        phase_delay = dm_out[0][1]/PSR_P0
        return DM_offset, phase_delay
    else:
        return 0, phase_off
    


def plot_OPT(save_path, archive_INIT, archive_OPT, fit_params):

    intensity_profile_INIT = archive_INIT.get_intensity_prof()
    intensity_profile_OPT = archive_OPT.get_intensity_prof()

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)

    prof_nbins = len(intensity_profile_OPT)
    phase = np.linspace(0, 1, prof_nbins)

    def func_r(x, A1, A2, x1, x2, d):
        g1 = profile_2((x-x1) % 1, *profile_pars)
        g2 = profile_2((x-x2) % 1, *profile_pars)
        return A1*g1 + A2*g2 + d

    archives = [archive_INIT, archive_OPT]
    titles = ['INIT', 'OPT']

    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(3, 2, hspace=0.0, wspace=0.25)

    axes = []
    first_ax = fig.add_subplot(gs[0, 0])
    axes.append([first_ax])

    for j in range(1, 2):
        ax = fig.add_subplot(gs[0, j], sharex=first_ax)
        axes[0].append(ax)

    for i in range(1, 3):
        row = []
        for j in range(2):
            ax = fig.add_subplot(gs[i, j], sharex=first_ax)
            row.append(ax)
        axes.append(row)

    for i in range(2):
        for j in range(2):
            axes[i][j].tick_params(labelbottom=False)

    for col, (archive, title) in enumerate(zip(archives, titles)):
        IP = archive.get_intensity_prof()
        FP = archive.get_freq_phase()
        TP = archive.get_time_phase()

        if np.max(FP) != 0:
            FP = FP / np.max(FP)

        phase = np.linspace(0, 1, len(IP))
        nchans = FP.shape[0]
        Tobs = TP.shape[0]

        axes[0][col].plot(phase, IP)
        axes[0][col].set_title(f'{title}, S/N: {archive.get_SNR():.2f}')
        axes[0][col].set_ylabel('Intensity')

        axes[1][col].imshow(
            FP,
            origin='lower',
            extent=[0, 1, 0, nchans],
            aspect='auto'
        )
        axes[1][col].set_ylabel('Channel number')

        axes[2][col].imshow(
            TP,
            origin='lower',
            extent=[0, 1, 0, Tobs],
            aspect='auto'
        )
        axes[2][col].set_ylabel('Time (s)')
        axes[2][col].set_xlabel('Phase')

    axes[0][1].plot(phase, np.roll(func_r(phase, *fit_params[0]), phase_corr-len(IP)//2), 'C1--')

    plt.savefig(save_path, dpi=200, bbox_inches="tight")

    

def plot_INIT(save_path, archive_INIT, out):
    fig = plt.figure(figsize=(12, 8))

    gs = fig.add_gridspec(3, 3, height_ratios=[1, 2, 2], hspace=0.0, wspace=0.3)

    axes = []
    first_ax = fig.add_subplot(gs[0, 0])
    axes.append([first_ax])

    for j in range(1, 3):
        ax = fig.add_subplot(gs[0, j], sharex=first_ax)
        axes[0].append(ax)

    for i in range(1, 3):
        row = []
        for j in range(3):
            ax = fig.add_subplot(gs[i, j], sharex=first_ax)
            row.append(ax)
        axes.append(row)

    for i in range(2):
        for j in range(3):
            axes[i][j].tick_params(labelbottom=False)

    freq_deriv, phase_offset, time, time_amp, Tobs = out

    IP = archive_INIT.get_intensity_prof()
    TP = archive_INIT.get_time_phase()
    FP = archive_INIT.get_freq_phase()

    nchans = len(FP)
    time_nbins = len(TP)
    phase_bins = len(IP)

    params, phase_corr, _ = get_IP_interp(IP)

    prof2D = scale_freq_phase(FP, IP)
    prof2D /= np.max(prof2D)
    prof2D = np.roll(prof2D, phase_bins//2+phase_corr, axis=1)
    FP /= np.max(FP)

    axes[1][0].imshow(FP, origin='lower', extent=[0, 1, 0, nchans], aspect='auto')
    axes[1][1].imshow(prof2D, origin='lower', extent=[0, 1, 0, nchans], aspect='auto')
    axes[1][2].imshow(FP-prof2D, origin='lower', extent=[0, 1, 0, nchans], aspect='auto')

    phase_arr = np.linspace(0, 1, phase_bins)
    axes[0][0].plot(phase_arr, IP)

    phase_plot = np.linspace(0, 1, phase_bins*10)
    axes[0][1].plot(phase_plot, np.roll(profile_2(phase_plot, *params), -10*(phase_bins//2-phase_corr)), 'C0-')
    axes[0][1].plot(phase_plot, np.roll(profile(phase_plot, *params[:3]), -10*(phase_bins//2-phase_corr)), 'C3--', lw=1)
    axes[0][1].plot(phase_plot, np.roll(profile(phase_plot, *params[3:]), -10*(phase_bins//2-phase_corr)), 'C3--', lw=1)

    axes[0][2].plot(phase_arr, IP-np.roll(profile_2(phase_arr, *params), phase_corr-phase_bins//2), 'C2-')

    for i in range(3):
        axes[0][i].set_ylabel('Intensity')
        axes[1][i].set_ylabel('Channel number')
        axes[2][i].set_ylabel('Time (s)')

    for i in range(3):
        axes[2][i].set_xlabel('Phase')
        axes[2][i].set_xlabel('Phase')
        axes[2][i].set_xlabel('Phase')

    axes[0][0].set_title('Pulsar fold')
    axes[0][1].set_title('Pulsar model')
    axes[0][2].set_title('Theoretical residuals')

    axes[2][0].imshow(TP, origin='lower', extent=[0, 1, 0, Tobs], aspect='auto')

    dt = Tobs / time_nbins
    time_tp = (np.arange(time_nbins) + 0.5) * dt - Tobs/2
    phase_tp = fit_f(time_tp, -phase_offset, *(-np.array([*freq_deriv.values()])))
    TP_model = np.zeros_like(TP)

    time_interp = PchipInterpolator(time, time_amp, extrapolate=True)
    for i in range(len(TP)):
        TP_model[i] = np.roll(profile_2(phase_arr-phase_tp[i], *params), -(phase_bins//2-phase_corr)) * time_interp(time_tp[i]+Tobs/2)

    axes[2][1].imshow(TP_model, origin='lower', extent=[0, 1, 0, Tobs], aspect='auto')
    axes[2][2].imshow(TP/np.mean(TP) - TP_model/np.mean(TP_model), origin='lower', extent=[0, 1, 0, Tobs], aspect='auto')

    plt.savefig(save_path, dpi=200, bbox_inches="tight")