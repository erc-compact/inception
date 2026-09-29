import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.interpolate import PchipInterpolator


DM_CONST_PULSARX = 1 / 2.41e-4


def bin_phase(nbins):
    return (np.arange(nbins) + 0.5) / nbins


def template_phase(nbins):
    return np.arange(nbins) / nbins


def wrap_phase(phase):
    return (phase + 0.5) % 1 - 0.5


def normalise_profile(prof):
    prof = prof - np.min(prof)
    peak = np.max(prof)
    return prof / peak if peak > 0 else prof


def normalise_light_curve(amps):
    amps = np.asarray(amps, dtype=np.float64)
    window = min(11, len(amps) - (1 - len(amps) % 2))
    if window > 3:
        amps = savgol_filter(amps, window_length=window, polyorder=3)
    amps = np.clip(amps, 0, None)
    return amps / np.mean(amps) if np.mean(amps) > 0 else np.ones_like(amps)


def profile(t, phase, sigma, Amp):
    g = norm(phase, sigma).pdf(t)
    return Amp * g

def profile_2(t, phase1, sigma1, Amp1, phase2, sigma2, Amp2):
    g1 = norm(phase1, sigma1).pdf(t)
    g2 = norm(phase2, sigma2).pdf(t)
    return Amp1 * g1 + Amp2 * g2


def two_pulse_model(profile_pars):
    def func_r(x, A1, A2, x1, x2, d):
        g1 = profile_2((x-x1) % 1, *profile_pars)
        g2 = profile_2((x-x2) % 1, *profile_pars)
        return A1*g1 + A2*g2 + d
    return func_r


def scale_freq_phase(freq_phase, intensity_profile):

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile)

    def fit_phase(phase, phase_off, Amp1, Amp2, base):
        g = profile_2((phase-phase_off)%1, *profile_pars[:2], Amp1, *profile_pars[3:5], Amp2)
        return g + base

    prof_nbins = freq_phase.shape[1]
    phase_arr = bin_phase(prof_nbins)
    phase_tmpl = template_phase(prof_nbins)
    freq_phase = np.roll(freq_phase, prof_nbins//2-phase_corr, axis=1)

    profile_s = []
    for time_phase_arr in freq_phase:
        if np.std(time_phase_arr) == 0:
            profile_s.append(np.zeros(prof_nbins))
            continue
        try:
            out = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 1, 0], bounds =[[-profile_pars[1], 0, 0, -np.inf], [profile_pars[1], np.inf, np.inf, np.inf]])
        except (RuntimeError, ValueError):
            profile_s.append(np.zeros(prof_nbins))
        else:
            updated_pars = profile_pars.copy()
            updated_pars[2] = out[0][1]
            updated_pars[5] = out[0][2]
            profile_s.append(profile_2((phase_tmpl-out[0][0]) % 1, *updated_pars))

    return np.array(profile_s)


def get_IP_interp(intensity_profile):

    true_prof = intensity_profile / np.max(intensity_profile)

    max_ind = np.argmax(true_prof)

    true_prof = np.roll(true_prof, len(true_prof)//2-max_ind)
    phase = bin_phase(len(true_prof))
    min_sigma = 0.25 / len(true_prof)
    bounds = [[0, min_sigma, 0, 0, min_sigma, 0, -1], [1, 1, 1, 1, 1, 1, 1]]

    def profile_2_base(t, phase1, sigma1, Amp1, phase2, sigma2, Amp2, base):
        return profile_2(t, phase1, sigma1, Amp1, phase2, sigma2, Amp2) + base

    starts = [[0.49, 0.05, 0.5, 0.51, 0.05, 0.5, 0]]
    for width in (0.015, 0.03, 0.06):
        for dx in (-0.06, -0.03, -0.01, 0.01, 0.03, 0.06):
            starts.append([0.5, width, 2.5*width, 0.5+dx, 1.5*width, 0.75*width, 0])

    out, best_cost = None, np.inf
    for p0 in starts:
        try:
            fit = curve_fit(profile_2_base, phase, true_prof, p0=p0, bounds=bounds)
        except (RuntimeError, ValueError):
            continue
        cost = np.sum((true_prof-profile_2_base(phase, *fit[0]))**2)
        if cost < best_cost:
            out, best_cost = fit, cost

    true_prof = true_prof - out[0][6]
    out = (out[0][:6], out[1][:6, :6])

    noise = true_prof-profile_2(phase, *out[0])
    SNR = np.sqrt(np.sum((true_prof-np.mean(noise))**2)/np.std(noise)**2)

    return out[0], max_ind, SNR


def fit_f(t, theta0, f0, f1):
    return theta0 + f0*t + (1/2)*f1*t**2


def fit_time_phase(time_phase, intensity_profile, obs_len, min_subints=5):

    profile_pars, phase_corr, SNR = get_IP_interp(intensity_profile)

    def fit_phase(phase, phase_off, Amp, base):
        g = profile_2((phase-phase_off)%1, *profile_pars[:2], profile_pars[2]*Amp, *profile_pars[3:5], profile_pars[5]*Amp)
        return g + base

    prof_nbins = len(intensity_profile)
    time_nbins = len(time_phase)
    phase_arr = bin_phase(prof_nbins)
    time_phase = np.roll(time_phase / np.max(time_phase), prof_nbins//2-phase_corr, axis=1)

    dt = obs_len / time_nbins
    time = (np.arange(time_nbins) + 0.5) * dt - obs_len/2

    good, phase, err = [], [], []
    for i, time_phase_arr in enumerate(time_phase):
        if np.std(time_phase_arr) == 0:
            continue
        try:
            popt, pcov = curve_fit(fit_phase, phase_arr, time_phase_arr, p0=[1e-6, 1, 0])
        except (RuntimeError, ValueError):
            continue
        perr = np.sqrt(np.diag(pcov))
        if np.all(np.isfinite(perr)) and (popt[1] > 3*perr[1]):
            good.append(i)
            phase.append(wrap_phase(popt[0]))
            err.append(perr[0])

    theta = np.zeros(3)
    if len(good) >= min_subints:
        try:
            theta = curve_fit(fit_f, time[good], phase, sigma=np.array(err), p0=[1e-3,1e-6,1e-10])[0]
        except (RuntimeError, ValueError):
            theta = np.zeros(3)

    phase_offset = -theta[0]
    phase_shift = (phase_corr - prof_nbins//2) / prof_nbins

    freq_deriv = {}
    for i, fx in enumerate(theta[1:]):
        freq_deriv[f'F{i}'] = -fx

    model_phase = fit_f(time, *theta)
    time_amp = np.zeros(time_nbins)
    for i, time_phase_arr in enumerate(time_phase):
        g = profile_2((phase_arr-model_phase[i]) % 1, *profile_pars)
        design = np.vstack([g, np.ones_like(g)]).T
        time_amp[i] = np.linalg.lstsq(design, time_phase_arr, rcond=None)[0][0]

    time_amp_smooth = normalise_light_curve(time_amp)

    return freq_deriv, phase_offset, phase_shift, SNR, time+obs_len/2, time_amp_smooth


def fit_phase_offset(intensity_profile_OPT, intensity_profile_INIT):

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)
    func_r = two_pulse_model(profile_pars)

    prof_nbins = len(intensity_profile_OPT)
    phase = bin_phase(prof_nbins)
    intensity_profile_OPT = normalise_profile(np.roll(intensity_profile_OPT, prof_nbins//2-phase_corr))

    p0 = [0.5, -0.5, phase[np.argmax(intensity_profile_OPT)]-0.5, phase[np.argmin(intensity_profile_OPT)]-0.5, 0.5]
    bounds = [[0, -np.inf, -0.5, -np.inf, -np.inf], [np.inf, 0, 0.5, np.inf, np.inf]]
    try:
        out = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0, bounds=bounds)
    except (RuntimeError, ValueError):
        return None

    phase_offset = wrap_phase(out[0][3]-out[0][2])
    SNR_scale = np.abs(out[0][0]/out[0][1])

    return phase_offset, SNR_scale, out


def fit_two_pulse(intensity_profile_OPT, func_r, phase_corr, x1, x2, width=0.15, SNR_limit=5):
    if np.std(intensity_profile_OPT) == 0:
        return None

    prof_nbins = len(intensity_profile_OPT)
    phase = bin_phase(prof_nbins)
    intensity_profile_OPT = normalise_profile(np.roll(intensity_profile_OPT, prof_nbins//2-phase_corr))

    p0 = [0.5, -0.5, x1, x2, np.clip(np.median(intensity_profile_OPT), 0.01, 0.99)]
    bounds = [[0, -1.5, x1-width, x2-width, -0.5], [1.5, 0, x1+width, x2+width, 1.5]]
    try:
        popt, pcov = curve_fit(func_r, phase, intensity_profile_OPT, p0=p0, bounds=bounds)
    except (RuntimeError, ValueError):
        return None

    residual_std = np.std(intensity_profile_OPT-func_r(phase, *popt))
    if (popt[0] < SNR_limit*residual_std) or (-popt[1] < SNR_limit*residual_std):
        return None

    offset_var = pcov[2, 2] + pcov[3, 3] - 2*pcov[2, 3]
    if (not np.isfinite(offset_var)) or (offset_var <= 0):
        return None

    return popt, np.sqrt(offset_var)


def fit_subint_phase_offset(time_phase_OPT, intensity_profile_INIT, fit_params):

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)
    func_r = two_pulse_model(profile_pars)
    x1, x2 = fit_params[0][2], fit_params[0][3]

    snr_arr = []
    for subint_i in time_phase_OPT:
        fit = fit_two_pulse(subint_i, func_r, phase_corr, x1, x2)
        snr_arr.append(np.abs(fit[0][0]/fit[0][1]) if fit else np.nan)

    snr_arr = np.array(snr_arr)
    if np.all(np.isnan(snr_arr)):
        return np.ones(len(snr_arr))
    snr_arr[np.isnan(snr_arr)] = np.nanmedian(snr_arr)

    return normalise_light_curve(snr_arr)


def fit_chan_phase_offset(freq_phase_OPT, intensity_profile_INIT, freq_arr, fit_params, PSR_P0, min_chans=8, sigma_level=3):

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)
    func_r = two_pulse_model(profile_pars)
    x1, x2 = fit_params[0][2], fit_params[0][3]
    phase_off = wrap_phase(x2-x1)
    f_top = np.max(freq_arr)

    def delay(freq_arr, DM, off):
        return -DM * DM_CONST_PULSARX * (1/freq_arr**2 - 1/f_top**2) + off

    freqs, delays, delay_err = [], [], []
    for freq, chan_i in zip(freq_arr, freq_phase_OPT):
        fit = fit_two_pulse(chan_i, func_r, phase_corr, x1, x2)
        if fit:
            freqs.append(freq)
            delays.append((fit[0][3]-fit[0][2]) * PSR_P0)
            delay_err.append(fit[1] * PSR_P0)

    if len(freqs) < max(min_chans, 3):
        return 0, phase_off

    try:
        dm_out = curve_fit(delay, np.array(freqs), np.array(delays), sigma=np.array(delay_err))
    except (RuntimeError, ValueError):
        return 0, phase_off

    DM_err = np.sqrt(dm_out[1][0, 0])
    if np.isfinite(DM_err) and (DM_err > 0) and (abs(dm_out[0][0]/DM_err) > sigma_level):
        return dm_out[0][0], wrap_phase(dm_out[0][1]/PSR_P0)
    else:
        return 0, phase_off



def plot_OPT(save_path, archive_INIT, archive_OPT, fit_params):

    intensity_profile_INIT = archive_INIT.get_intensity_prof()
    intensity_profile_OPT = archive_OPT.get_intensity_prof()

    profile_pars, phase_corr, _ = get_IP_interp(intensity_profile_INIT)
    func_r = two_pulse_model(profile_pars)

    prof_nbins = len(intensity_profile_OPT)
    phase = bin_phase(prof_nbins)

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
        IP = normalise_profile(archive.get_intensity_prof())
        FP = archive.get_freq_phase()
        TP = archive.get_time_phase()

        if np.max(FP) != 0:
            FP = FP / np.max(FP)

        phase = bin_phase(len(IP))
        nchans = FP.shape[0]
        Tobs = TP.shape[0]

        axes[0][col].plot(phase, IP)
        axes[0][col].set_title(f'{title}, S/N: {archive.get_SNR():.2f}')
        axes[0][col].set_ylabel('Intensity')

        axes[1][col].imshow(
            FP,
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
    plt.close(fig)



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
    IP = IP / np.max(IP)
    TP = archive_INIT.get_time_phase()
    FP = archive_INIT.get_freq_phase()

    nchans = len(FP)
    time_nbins = len(TP)
    phase_bins = len(IP)

    params, phase_corr, _ = get_IP_interp(IP)

    prof2D = scale_freq_phase(FP, IP)
    prof2D = 0.5 * (prof2D + np.roll(prof2D, -1, axis=1))
    prof2D /= np.max(prof2D)
    prof2D = np.roll(prof2D, phase_bins//2+phase_corr, axis=1)
    FP /= np.max(FP)

    axes[1][0].imshow(FP,  extent=[0, 1, 0, nchans], aspect='auto')
    axes[1][1].imshow(prof2D,  extent=[0, 1, 0, nchans], aspect='auto')
    axes[1][2].imshow(FP-prof2D,  extent=[0, 1, 0, nchans], aspect='auto')

    phase_arr = bin_phase(phase_bins)
    axes[0][0].plot(phase_arr, IP)

    phase_plot = (np.arange(phase_bins*10) + 0.5) / (phase_bins*10)
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
        TP_model[i] = np.roll(profile_2((phase_arr-phase_tp[i]) % 1, *params), -(phase_bins//2-phase_corr)) * time_interp(time_tp[i]+Tobs/2)

    axes[2][1].imshow(TP_model, origin='lower', extent=[0, 1, 0, Tobs], aspect='auto')
    axes[2][2].imshow(TP/np.mean(TP) - TP_model/np.mean(TP_model), origin='lower', extent=[0, 1, 0, Tobs], aspect='auto')

    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
