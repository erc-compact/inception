import os
import re
import sys
import numpy as np
import astropy.units as u
from math import factorial
from decimal import Decimal
from functools import lru_cache
from astropy.time import TimeDelta
from scipy.interpolate import interp1d

from .propagation_effects import TEMPO_DM_CONST


@lru_cache(maxsize=None)
def load_pint():
    try:
        import pint.logging as logging      # type: ignore
        logging.setup('ERROR')
        import pint.models as models        # type: ignore
        from pint.polycos import Polycos    # type: ignore
    except ImportError:
        sys.exit('pint-pulsar package not installed, cannot use pint_polycos.')
    return models, Polycos


def read_par(par_path):
    with open(par_path) as par_file:
        return par_file.read().splitlines()


def get_par_value(lines, key):
    for line in lines:
        parts = line.split()
        if len(parts) > 1 and parts[0] == key:
            return Decimal(parts[1].replace('D', 'E').replace('d', 'e'))
    return None


def set_par_value(lines, key, value):
    for n, line in enumerate(lines):
        parts = line.split()
        if len(parts) > 1 and parts[0] == key:
            if Decimal(parts[1].replace('D', 'E').replace('d', 'e')) == value:
                return
            token = format(value, 'E') if re.search('[EeDd]', parts[1]) else format(value, 'f')
            head, tail = re.match(r'(\s*\S+\s+)\S+(.*)', line).groups()
            lines[n] = head + token + tail
            return
    if value:
        lines.append(f'{key:<15}{value}')


def par_frequency(lines):
    f0 = get_par_value(lines, 'F0')
    if f0:
        return float(f0)
    p0 = get_par_value(lines, 'P0')
    return 1/float(p0) if p0 else 0.


def ephemeris_values(ephemeris_path):
    if (not ephemeris_path) or (not os.path.exists(ephemeris_path)):
        return 0., 0.
    if ephemeris_path.endswith('.par'):
        lines = read_par(ephemeris_path)
        return par_frequency(lines), float(get_par_value(lines, 'DM') or 0)
    _, Polycos = load_pint()
    polyco_table = Polycos.read(ephemeris_path).polycoTable
    return float(polyco_table['entry'][0].f0), float(polyco_table['dm'][0])


def solar_system_ephemeris(ephem):
    if not os.path.isfile(ephem):
        return ephem
    from pint.solar_system_ephemerides import load_kernel # type: ignore
    ephem_name = os.path.splitext(os.path.basename(ephem))[0].lower()
    load_kernel(ephem_name, path=ephem)
    return ephem_name


def create_polycos(par_file, pulsar_pars, obs, output_path):
    models, Polycos = load_pint()
    timing_model = models.get_model(par_file, EPHEM=solar_system_ephemeris(obs.ephem))

    t_mid = obs.obs_start + obs.obs_len/2 * u.s.to(u.day)
    polco_range = obs.obs_len/2 + 100*u.min.to(u.s)
    start, end = t_mid - polco_range*u.s.to(u.day), t_mid + polco_range*u.s.to(u.day)

    polycos_coeff = max(1, pulsar_pars['pint_N'])
    polycos_tspan = max(1, pulsar_pars['pint_T']) # minutes
    gen_poly = Polycos.generate_polycos(timing_model, start, end, obs.tempo_id,
                                        polycos_tspan, polycos_coeff,
                                        obs.f0, progress=False)

    polycos_path = f"{output_path}/{pulsar_pars['ID']}.polycos"
    gen_poly.write_polyco_file(polycos_path)
    return polycos_path


def write_parfile(pulsar_model, par_path, new_path):
    lines = read_par(par_path)

    prop = pulsar_model.prop_effect
    DM_offset = Decimal(repr(float(pulsar_model.pulsar_pars['DM'] * prop.DM_const / TEMPO_DM_CONST)))
    set_par_value(lines, 'DM', (get_par_value(lines, 'DM') or Decimal(0)) + DM_offset)

    epoch = pulsar_model.pepoch + (pulsar_model.accepoch if pulsar_model.AX_list else 0) * u.s.to(u.day)
    if pulsar_model.pulsar_pars['DM_ref'] != 'inf':
        epoch -= prop.DM * prop.DM_const / pulsar_model.obs.high_f**2 * u.s.to(u.day)
    if pulsar_model.pulsar_pars['frame'] == 'topo':
        epoch = float(pulsar_model.obs.topo2bary([epoch])[0])

    par_pepoch = get_par_value(lines, 'PEPOCH')
    if par_pepoch is None:
        par_pepoch = Decimal(repr(float(epoch)))
        set_par_value(lines, 'PEPOCH', par_pepoch)
    shift = (float(par_pepoch) - epoch) * u.day.to(u.s)

    FX = pulsar_model.FX_doppler
    deltas = [sum(FX[k+j] * shift**j / factorial(j) for j in range(len(FX) - k)) for k in range(len(FX))]
    highest = max([k for k, delta in enumerate(deltas) if delta], default=-1)
    for k in range(highest + 1):
        value = get_par_value(lines, f'F{k}')
        set_par_value(lines, f'F{k}', (value or Decimal(0)) + Decimal(repr(float(deltas[k]))))

    with open(new_path, 'w') as par_file:
        par_file.write('\n'.join(lines) + '\n')
    return new_path


class PolycoPhase:
    def __init__(self, polycos_path, obs, prop_effect, DM_ref):
        _, Polycos = load_pint()
        polyco_table = Polycos.read(polycos_path).polycoTable

        self.ref_freq = float(polyco_table['obsfreq'][0])
        self.ephemeris_DM = float(polyco_table['dm'][0])
        self.inv_ref = 0. if DM_ref == 'inf' else 1/obs.high_f**2
        self.DM_sweep = prop_effect.DM * prop_effect.DM_const
        self.chan_delays = self.delay(obs.freq_arr)

        cover_lo = (np.min(np.asarray(polyco_table['t_start'], dtype=np.float64)) - obs.obs_start) * 86400
        cover_hi = (np.max(np.asarray(polyco_table['t_stop'], dtype=np.float64)) - obs.obs_start) * 86400
        interp_lo = max(cover_lo, min(0., np.min(self.chan_delays)) - 3600)
        interp_hi = min(cover_hi, obs.obs_len + obs.dt + max(0., np.max(self.chan_delays)) + 3600)
        interp_sec = np.linspace(interp_lo, interp_hi, 10**5)
        interp_time = obs.obs_start_time_TIME + TimeDelta(interp_sec, format='sec')
        interp_mjd = (interp_time.jd1 - 2400000.5).astype(np.longdouble) + interp_time.jd2.astype(np.longdouble)
        pulse_int, pulse_frac = self.entry_phase(polyco_table, interp_mjd)
        rel_phase = (pulse_int - pulse_int[0]).astype(np.float64) + pulse_frac.astype(np.float64)

        self.interp = interp1d(interp_sec, rel_phase)

    @staticmethod
    def entry_phase(polyco_table, mjd):
        tmids = np.asarray(polyco_table['tmid'], dtype=np.longdouble)
        order = np.argsort(tmids)
        nearest = order[np.searchsorted((tmids[order][1:] + tmids[order][:-1]) / 2, mjd)]
        pulse_int = np.zeros(len(mjd), dtype=np.longdouble)
        pulse_frac = np.zeros(len(mjd), dtype=np.longdouble)
        for k in np.unique(nearest):
            in_entry = nearest == k
            phase = polyco_table['entry'][k].evalabsphase(mjd[in_entry])
            pulse_int[in_entry] = phase.int.value
            pulse_frac[in_entry] = phase.frac.value
        return pulse_int, pulse_frac

    def delay(self, freq):
        ephemeris_shift = self.ephemeris_DM * TEMPO_DM_CONST * (1/self.ref_freq**2 - self.inv_ref)
        return ephemeris_shift - self.DM_sweep * (1/np.asarray(freq, dtype=np.float64)**2 - self.inv_ref)

    def __call__(self, t_sec, freq):
        return self.interp(t_sec + self.delay(freq))

    def channels(self, t_sec):
        return self.interp(t_sec[:, None] + self.chan_delays[None, :])
