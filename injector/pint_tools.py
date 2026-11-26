import os
import re
import sys
import numpy as np
import astropy.units as u
from math import comb, factorial
from decimal import Decimal
from functools import lru_cache
from astropy.time import TimeDelta
from scipy.interpolate import interp1d, CubicSpline

from .propagation_effects import TEMPO_DM_CONST

MIN_PER_DAY = 1440
POLYCO_TOLERANCE = 1e-6
DELAY_GRID_STEP = 600.


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
    if not ephemeris_path:
        return 0., 0.
    if not os.path.exists(ephemeris_path):
        sys.exit(f'Ephemeris file not found: {ephemeris_path}')
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


def shift_polyco_phases(tables):
    from pint.phase import Phase    # type: ignore
    entries = [entry for table in tables for entry in table['entry']]
    lowest = min(int(entry.rphase.int[0]) for entry in entries)
    if lowest < 1:
        shift = Phase(np.longdouble(1 - lowest))
        for entry in entries:
            entry.rphase = entry.rphase + shift


def anchor_written_tmid(table):
    from pint.phase import Phase    # type: ignore
    for row in range(len(table)):
        entry = table['entry'][row]
        written = np.longdouble('{:20.11f}'.format(table['tmid'][row]).strip())
        delta = (written - entry.tmid.value) * MIN_PER_DAY
        coeffs = entry.coeffs
        n_coeff = len(coeffs)
        shifted = np.array([sum(comb(i, k) * coeffs[i] * delta**(i - k) for i in range(k, n_coeff))
                            for k in range(n_coeff)], dtype=np.longdouble)
        entry.rphase = entry.rphase + Phase(shifted[0] - coeffs[0] + 60 * entry.f0 * delta)
        shifted[0] = coeffs[0]
        entry.coeffs = shifted
        entry.tmid = written * u.day
        entry.tstart = entry.tmid - entry.mjdspan / 2
        entry.tstop = entry.tmid + entry.mjdspan / 2
        table['tmid'][row] = written


def polyco_frequencies(obs, n_freq):
    centre = round(float(obs.f0), 3)
    if n_freq <= 1:
        return np.array([centre])
    inv_sq = np.linspace(1/obs.high_f**2, 1/obs.low_f**2, n_freq)
    return np.unique(np.round(np.append(1/np.sqrt(inv_sq), centre), 3))


def generate_polyco_table(timing_model, Polycos, obs, freq, n_coeff, t_span, margin):
    t_mid = obs.obs_start + obs.obs_len/2 * u.s.to(u.day)
    half = (obs.obs_len/2 + margin) * u.s.to(u.day)
    polycos = Polycos.generate_polycos(timing_model, t_mid - half, t_mid + half, obs.tempo_id,
                                       t_span, n_coeff, freq, progress=False)
    return polycos.polycoTable


def write_polyco_tables(tables, path):
    from pint.polycos import tempo_polyco_table_writer    # type: ignore
    with open(path, 'w') as file:
        for k, table in enumerate(tables):
            part = f'{path}.part{k}'
            tempo_polyco_table_writer(table, part)
            with open(part) as part_file:
                file.write(part_file.read())
            os.remove(part)


def polyco_error(timing_model, Polycos, path, obs, freqs, n_test=41):
    import pint.toa as toa    # type: ignore
    table = Polycos.read(path).polycoTable
    obsfreq = np.asarray(table['obsfreq'], dtype=np.float64)

    times = obs.obs_start_time_TIME + TimeDelta(np.linspace(0, obs.obs_len, n_test), format='sec')
    mjd = (times.jd1 - 2400000.5).astype(np.longdouble) + times.jd2.astype(np.longdouble)
    worst = 0.
    for freq in freqs:
        toas = toa.get_TOAs_array(times, obs=obs.tempo_id, freqs=freq * u.MHz, ephem=timing_model.EPHEM.value)
        model_frac = timing_model.phase(toas, abs_phase=True).frac.value.astype(np.float64)
        _, polyco_frac = PolycoPhase.entry_phase(table[obsfreq == freq], mjd)
        diff = polyco_frac.astype(np.float64) - model_frac
        worst = max(worst, float(np.max(np.abs(diff - np.round(diff)))))
    return worst


def create_polycos(par_file, pulsar_pars, obs, output_path):
    models, Polycos = load_pint()
    timing_model = models.get_model(par_file, EPHEM=solar_system_ephemeris(obs.ephem))
    if 'AbsPhase' not in timing_model.components:
        print(f'WARNING: {par_file} has no TZRMJD. PINT sets its own phase reference for this observation, '
              'so injections into different epochs will not be phase-connected.', flush=True)

    freqs = polyco_frequencies(obs, pulsar_pars.get('pint_nfreq', 8))
    n_coeff = max(1, pulsar_pars['pint_N'])
    t_span = max(1, int(pulsar_pars['pint_T']))
    polycos_path = f"{output_path}/{pulsar_pars['ID']}.polycos"

    while True:
        tables = [generate_polyco_table(timing_model, Polycos, obs, freq, n_coeff, t_span, 100*u.min.to(u.s))
                  for freq in freqs]
        for table in tables:
            anchor_written_tmid(table)
        shift_polyco_phases(tables)
        write_polyco_tables(tables, polycos_path)
        error = polyco_error(timing_model, Polycos, polycos_path, obs, freqs)
        if error <= POLYCO_TOLERANCE or t_span == 1:
            break
        t_span = max(1, t_span // 2)

    print(f"PINT polycos for {pulsar_pars['ID']}: {len(freqs)} frequencies, {t_span}-min segments, "
          f"largest error {error:.1e} turns", flush=True)
    if error > POLYCO_TOLERANCE:
        print(f"WARNING: PINT polycos for {pulsar_pars['ID']} deviate from the timing model by up to "
              f"{error:.1e} turns; increase pint_N.", flush=True)
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
        obsfreq = np.asarray(polyco_table['obsfreq'], dtype=np.float64)
        nodes = np.unique(obsfreq)

        self.obs = obs
        self.ref_freq = nodes[np.argmin(np.abs(nodes - obs.f0))]
        self.ref_table = polyco_table[obsfreq == self.ref_freq]
        self.ephemeris_DM = float(polyco_table['dm'][0])
        self.inv_ref = 0. if DM_ref == 'inf' else 1/obs.high_f**2
        self.offset_sweep = prop_effect.DM * prop_effect.DM_const - self.ephemeris_DM * TEMPO_DM_CONST

        n_grid = max(2, int(np.ceil(obs.obs_len / DELAY_GRID_STEP)) + 1)
        self.t_grid = np.linspace(0., obs.obs_len, n_grid)
        self.node_delays = self.measure_node_delays(polyco_table, obsfreq, nodes) if len(nodes) > 1 else None
        self.chan_table = self.delay_table(obs.freq_arr)

        cover_lo = (np.min(np.asarray(self.ref_table['t_start'], dtype=np.float64)) - obs.obs_start) * 86400
        cover_hi = (np.max(np.asarray(self.ref_table['t_stop'], dtype=np.float64)) - obs.obs_start) * 86400
        interp_lo = max(cover_lo, min(0., np.min(self.chan_table)) - 3600)
        interp_hi = min(cover_hi, obs.obs_len + obs.dt + max(0., np.max(self.chan_table)) + 3600)
        interp_sec = np.linspace(interp_lo, interp_hi, 10**5)
        pulse_int, pulse_frac = self.entry_phase(self.ref_table, self.to_mjd(interp_sec))
        rel_phase = (pulse_int - pulse_int[0]).astype(np.float64) + pulse_frac.astype(np.float64)

        self.interp = interp1d(interp_sec, rel_phase)

    def to_mjd(self, t_sec):
        times = self.obs.obs_start_time_TIME + TimeDelta(t_sec, format='sec')
        return (times.jd1 - 2400000.5).astype(np.longdouble) + times.jd2.astype(np.longdouble)

    @staticmethod
    def nearest_entries(polyco_table, mjd):
        tmids = np.asarray(polyco_table['tmid'], dtype=np.longdouble)
        order = np.argsort(tmids)
        return order[np.searchsorted((tmids[order][1:] + tmids[order][:-1]) / 2, mjd)]

    @staticmethod
    def entry_phase(polyco_table, mjd):
        nearest = PolycoPhase.nearest_entries(polyco_table, mjd)
        pulse_int = np.zeros(len(mjd), dtype=np.longdouble)
        pulse_frac = np.zeros(len(mjd), dtype=np.longdouble)
        for k in np.unique(nearest):
            in_entry = nearest == k
            phase = polyco_table['entry'][k].evalabsphase(mjd[in_entry])
            pulse_int[in_entry] = phase.int.value
            pulse_frac[in_entry] = phase.frac.value
        return pulse_int, pulse_frac

    @staticmethod
    def entry_freq(polyco_table, mjd):
        nearest = PolycoPhase.nearest_entries(polyco_table, mjd)
        freq = np.zeros(len(mjd), dtype=np.float64)
        for k in np.unique(nearest):
            in_entry = nearest == k
            freq[in_entry] = np.asarray(polyco_table['entry'][k].evalfreq(mjd[in_entry]), dtype=np.float64)
        return freq

    def measure_node_delays(self, polyco_table, obsfreq, nodes):
        mjd = self.to_mjd(self.t_grid)
        delays = np.zeros((len(nodes), len(self.t_grid)))
        for k, freq in enumerate(nodes):
            if freq == self.ref_freq:
                continue
            node_int, node_frac = self.entry_phase(polyco_table[obsfreq == freq], mjd)
            for _ in range(4):
                shifted = mjd + (delays[k] / 86400).astype(np.longdouble)
                ref_int, ref_frac = self.entry_phase(self.ref_table, shifted)
                residual = (node_int - ref_int).astype(np.float64) + (node_frac - ref_frac).astype(np.float64)
                delays[k] += residual / self.entry_freq(self.ref_table, shifted)
        remainder = delays - self.dispersion(nodes)[:, None]
        return CubicSpline(np.log(nodes), remainder, axis=0)

    def dispersion(self, freq):
        return self.ephemeris_DM * TEMPO_DM_CONST * (1/self.ref_freq**2 - 1/freq**2)

    def delay_table(self, freq):
        freq = np.asarray(freq, dtype=np.float64)
        chromatic = np.broadcast_to(self.dispersion(freq), (len(self.t_grid),) + freq.shape)
        if self.node_delays is not None:
            chromatic = chromatic + np.moveaxis(self.node_delays(np.log(freq)), -1, 0)
        return chromatic - self.offset_sweep * (1/freq**2 - self.inv_ref)

    def grid_weights(self, t_sec):
        position = np.interp(t_sec, self.t_grid, np.arange(len(self.t_grid)))
        lower = np.minimum(position.astype(np.int64), len(self.t_grid) - 2)
        return lower, position - lower

    def __call__(self, t_sec, freq):
        t_sec = np.asarray(t_sec, dtype=np.float64)
        t_flat, f_flat = (np.ravel(a) for a in np.broadcast_arrays(t_sec, freq))
        table = self.delay_table(f_flat)
        lower, weight = self.grid_weights(t_flat)
        columns = np.arange(len(f_flat))
        delay = table[lower, columns] * (1 - weight) + table[lower + 1, columns] * weight
        return self.interp(t_flat + delay).reshape(np.broadcast(t_sec, freq).shape)

    def channels(self, t_sec):
        lower, weight = self.grid_weights(t_sec)
        delay = self.chan_table[lower] * (1 - weight)[:, None] + self.chan_table[lower + 1] * weight[:, None]
        return self.interp(t_sec[:, None] + delay)
