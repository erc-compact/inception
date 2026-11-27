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
MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
POLYCO_TOLERANCE = 1e-6
BAND_TOLERANCE = 1e-8
MAX_NFREQ = 64
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


def par_ephemeris(lines):
    for line in lines:
        parts = line.split()
        if len(parts) > 1 and parts[0] == 'EPHEM':
            return parts[1].lower()
    return 'de421'


def set_par_ephemeris(lines, name):
    for n, line in enumerate(lines):
        parts = line.split()
        if len(parts) > 1 and parts[0] == 'EPHEM':
            head, tail = re.match(r'(\s*\S+\s+)\S+(.*)', line).groups()
            lines[n] = head + name.upper() + tail
            return
    lines.append(f'{"EPHEM":<15}{name.upper()}')


@lru_cache(maxsize=None)
def pint_ephemeris(ephem, par_ephem):
    # PINT follows the par file's EPHEM; the --ephem file supplies it when the names match. Returns the kernel
    # name for PINT and the value astropy needs while PINT runs (PINT reads astropy's global ephemeris).
    from astropy.coordinates import solar_system_ephemeris    # type: ignore
    from pint.solar_system_ephemerides import load_kernel    # type: ignore
    previous = solar_system_ephemeris.get()
    given = os.path.splitext(os.path.basename(ephem))[0].lower() if os.path.isfile(ephem) else None

    def load(name, path=None):
        loaded = load_kernel(name, path=path)
        return loaded if isinstance(loaded, (str, os.PathLike)) else solar_system_ephemeris.get()

    try:
        if given == par_ephem:
            return par_ephem, load(par_ephem, ephem)
        try:
            return par_ephem, load(par_ephem)
        except Exception:
            pass
        if given:
            print(f'WARNING: PINT could not load {par_ephem.upper()} named in the par file; using the --ephem kernel '
                  f'{given.upper()} instead (written to the output par file).', flush=True)
            return given, load(given, ephem)
        print(f'WARNING: PINT could not load {par_ephem.upper()} named in the par file and --ephem gives no kernel; '
              "using astropy's built-in ephemeris, good to about 10 us only (not timing-grade). Give --ephem "
              f'with the {par_ephem} kernel.', flush=True)
        return 'builtin', 'builtin'
    finally:
        solar_system_ephemeris.set(previous)


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


BIPM_FALLBACK = set()


def clock_settings(timing_model):
    # PINT's get_TOAs reads CLOCK from the model, get_TOAs_array (and generate_polycos) does not
    clock = getattr(timing_model, 'CLOCK', None)
    clock = str(clock.value if clock is not None and clock.value is not None else '').replace(' ', '').upper()
    if clock in ('TT(TAI)', 'UNCORR'):
        return {'include_bipm': False}
    if clock.startswith('TT(BIPM') and clock != 'TT(BIPM)' and clock[3:-1] not in BIPM_FALLBACK:
        return {'include_bipm': True, 'bipm_version': clock[3:-1]}
    return {'include_bipm': True}


def clock_label(timing_model):
    clock = getattr(timing_model, 'CLOCK', None)
    label = clock.value if clock is not None and clock.value else 'default clock'
    if 'bipm_version' not in clock_settings(timing_model) and str(label).upper().replace(' ', '')[3:-1] in BIPM_FALLBACK:
        label += ' unavailable, PINT default used'
    return label


def planet_shapiro(timing_model):
    param = getattr(timing_model, 'PLANET_SHAPIRO', None)
    return bool(param is not None and param.value)


def model_toas(timing_model, mjd, obs_code, freq):
    import pint.toa as toa    # type: ignore
    day = np.floor(mjd)
    times = (day.astype(np.float64), (mjd - day).astype(np.float64))
    common = dict(obs=obs_code, freqs=freq * u.MHz, ephem=timing_model.EPHEM.value, planets=planet_shapiro(timing_model))
    settings = clock_settings(timing_model)
    try:
        return toa.get_TOAs_array(times, **common, **settings)
    except Exception as error:
        version = settings.get('bipm_version')
        if not version or 'bipm' not in str(error).lower():
            raise
    # the TT(BIPMxxxx) realisations differ by about a nanosecond; a missing clock file should not stop the run
    BIPM_FALLBACK.add(version)
    print(f'WARNING: PINT cannot load the clock file for TT({version}) named in the par file; using PINT\'s default '
          f'TT(BIPM) realisation instead (they differ by about a nanosecond). Put tai2tt_{version.lower()}.clk into '
          'the cache for an exact match.', flush=True)
    return toa.get_TOAs_array(times, **common, **clock_settings(timing_model))


def model_phase(timing_model, mjd, obs_code, freq):
    mjd = np.asarray(mjd, dtype=np.longdouble)
    toas = model_toas(timing_model, mjd, obs_code, freq)
    phase = timing_model.phase(toas, abs_phase=True)
    order = np.argsort(np.asarray(toas.table['index']))
    return (np.asarray(phase.int.value, dtype=np.longdouble)[order],
            np.asarray(phase.frac.value, dtype=np.longdouble)[order], toas, order)


def coverage_margin(timing_model, pulsar_pars, obs):
    DM = abs(float(timing_model.DM.value or 0)) + abs(float(pulsar_pars.get('DM') or 0))
    return 600. + 1.5 * DM * TEMPO_DM_CONST / obs.low_f**2


def polyco_tmids(obs, t_span, margin):
    start = np.floor((obs.obs_start - margin / 86400) * MIN_PER_DAY)
    stop = np.ceil((obs.obs_start + (obs.obs_len + margin) / 86400) * MIN_PER_DAY)
    n_seg = max(1, int(np.ceil((stop - start) / t_span)))
    return (np.longdouble(start) + t_span * (np.arange(n_seg) + 0.5)) / MIN_PER_DAY


def orbit_columns(timing_model, toas, order, n_seg):
    if not timing_model.is_binary:
        return None
    try:
        orbital = np.asarray(timing_model.orbital_phase(toas, radians=False), dtype=np.float64)[order][:n_seg]
        component = timing_model.get_components_by_category()['pulsar_system'][0]
        PB = getattr(component, 'PB', None)
        f_orbit = 1/float(PB.value) if PB is not None and PB.value else float(component.FB0.value) * 86400
        return orbital, f_orbit
    except Exception:
        return None


def generate_polyco_table(timing_model, obs, freq, n_coeff, t_span, margin, n_nodes=20):
    # TEMPO-style polycos as PINT's generate_polycos makes them, but with one TOA call for all segments,
    # the par file's CLOCK and PLANET_SHAPIRO, and segments that cover the observation plus margin
    from astropy.table import Table
    from astropy.time import Time
    from pint.polycos import PolycoEntry    # type: ignore
    n_nodes = max(n_nodes, n_coeff + 1)
    tmids = polyco_tmids(obs, t_span, margin)
    n_seg = len(tmids)
    dt = np.linspace(-t_span/2, t_span/2, n_nodes)
    nodes = tmids[:, None] + dt.astype(np.longdouble)[None, :] / MIN_PER_DAY

    ints, fracs, toas, order = model_phase(timing_model, np.concatenate([tmids, nodes.ravel()]), obs.tempo_id, freq)
    ref_int, ref_frac = ints[:n_seg], fracs[:n_seg]
    f0 = timing_model.F0.value
    reduced = ((ints[n_seg:].reshape(n_seg, n_nodes) - ref_int[:, None]) - dt[None, :] * 60 * f0
               + (fracs[n_seg:].reshape(n_seg, n_nodes) - ref_frac[:, None]))
    coeffs = np.polyfit(dt, reduced.astype(np.float64).T, n_coeff - 1)[::-1].T
    orbit = orbit_columns(timing_model, toas, order, n_seg)

    rows = []
    for k, iso in enumerate(Time(tmids.astype(np.float64), format='mjd', scale='utc').iso):
        date, hms = iso.split()
        yy, mm, dd = date.split('-')
        entry = PolycoEntry(tmids[k], t_span, np.array([ref_int[k]]), np.array([ref_frac[k]]), f0, n_coeff, coeffs[k])
        row = {'psr': timing_model.PSR.value, 'date': f'{dd}-{MONTHS[int(mm) - 1]}-{yy[-2:]}',
               'utc': float(hms.replace(':', '')), 'tmid': tmids[k], 'dm': timing_model.DM.value,
               'doppler': 0., 'logrms': 0., 'mjd_span': t_span, 't_start': entry.tstart, 't_stop': entry.tstop,
               'obs': obs.tempo_id, 'obsfreq': float(freq)}
        if orbit is not None:
            row['binary_phase'], row['f_orbit'] = orbit[0][k], orbit[1]
        row['entry'] = entry
        rows.append(row)
    return Table(rows=rows, meta={'name': 'Polyco Data Table'})


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
    table = Polycos.read(path).polycoTable
    obsfreq = np.asarray(table['obsfreq'], dtype=np.float64)

    times = obs.obs_start_time_TIME + TimeDelta(np.linspace(0, obs.obs_len, n_test), format='sec')
    mjd = (times.jd1 - 2400000.5).astype(np.longdouble) + times.jd2.astype(np.longdouble)
    worst = 0.
    for freq in freqs:
        model_frac = model_phase(timing_model, mjd, obs.tempo_id, freq)[1].astype(np.float64)
        _, polyco_frac = PolycoPhase.entry_phase(table[obsfreq == freq], mjd)
        diff = polyco_frac.astype(np.float64) - model_frac
        worst = max(worst, float(np.max(np.abs(diff - np.round(diff)))))
    return worst


def band_error(timing_model, path, obs, n_times=9, max_chans=64):
    # injected phase (polycos + interpolation between their frequencies) against PINT, over the channels
    from types import SimpleNamespace
    n_chan = len(obs.freq_arr)
    chans = np.unique(np.round(np.linspace(0, n_chan - 1, min(max_chans, n_chan))).astype(int))
    freqs = np.asarray(obs.freq_arr, dtype=np.float64)[chans]
    t_sec = np.repeat(np.linspace(0., obs.obs_len, n_times), len(freqs))
    f_all = np.tile(freqs, n_times)

    prop = SimpleNamespace(DM=float(timing_model.DM.value or 0), DM_const=TEMPO_DM_CONST)
    injected = PolycoPhase(path, obs, prop, 'inf')(t_sec, f_all).astype(np.longdouble)
    times = obs.obs_start_time_TIME + TimeDelta(t_sec, format='sec')
    mjd = (times.jd1 - 2400000.5).astype(np.longdouble) + times.jd2.astype(np.longdouble)
    model_int, model_frac, _, _ = model_phase(timing_model, mjd, obs.tempo_id, f_all)

    diff = (model_int - injected + model_frac).astype(np.float64)
    seconds = (diff - np.round(diff)) / float(timing_model.F0.value)
    return len(freqs), float(np.max(np.abs(seconds))), float(np.sqrt(np.mean(seconds**2)))


def create_polycos(par_file, pulsar_pars, obs, output_path):
    from astropy.coordinates import solar_system_ephemeris    # type: ignore
    models, Polycos = load_pint()
    ephem_name, kernel = pint_ephemeris(str(obs.ephem), par_ephemeris(read_par(par_file)))
    with solar_system_ephemeris.set(kernel):
        return make_polycos(models, Polycos, par_file, pulsar_pars, obs, output_path, ephem_name)


def make_polycos(models, Polycos, par_file, pulsar_pars, obs, output_path, ephem_name):
    timing_model = models.get_model(par_file, EPHEM=ephem_name)
    if 'AbsPhase' not in timing_model.components:
        print(f'WARNING: {par_file} has no TZRMJD. PINT sets its own phase reference for this observation, '
              'so injections into different epochs will not be phase-connected.', flush=True)

    n_freq = max(1, int(pulsar_pars.get('pint_nfreq', 8)))
    n_coeff = max(1, pulsar_pars['pint_N'])
    t_span = max(1, int(pulsar_pars['pint_T']))
    polycos_path = f"{output_path}/{pulsar_pars['ID']}.polycos"
    margin = coverage_margin(timing_model, pulsar_pars, obs)

    while True:
        freqs = polyco_frequencies(obs, n_freq)
        while True:
            tables = [generate_polyco_table(timing_model, obs, freq, n_coeff, t_span, margin) for freq in freqs]
            for table in tables:
                anchor_written_tmid(table)
            shift_polyco_phases(tables)
            write_polyco_tables(tables, polycos_path)
            error = polyco_error(timing_model, Polycos, polycos_path, obs, freqs)
            if error <= POLYCO_TOLERANCE or t_span == 1:
                break
            t_span = max(1, t_span // 2)

        # pint_nfreq: 1 keeps the band centre only on purpose; otherwise add frequencies until the band is followed
        n_checked, largest, rms = band_error(timing_model, polycos_path, obs)
        if largest <= BAND_TOLERANCE or n_freq == 1 or n_freq >= MAX_NFREQ:
            break
        n_freq = min(MAX_NFREQ, 2 * n_freq)

    print(f"PINT polycos for {pulsar_pars['ID']} ({ephem_name.upper()}, {clock_label(timing_model)}): {len(freqs)} frequencies, "
          f"{t_span}-min segments. Injected phase against PINT over {n_checked} of {len(obs.freq_arr)} channels: "
          f"largest {largest*1e9:.1f} ns, rms {rms*1e9:.1f} ns", flush=True)
    if error > POLYCO_TOLERANCE:
        print(f"WARNING: PINT polycos for {pulsar_pars['ID']} deviate from the timing model by up to "
              f"{error:.1e} turns; increase pint_N.", flush=True)
    if largest > BAND_TOLERANCE:
        hint = 'pint_nfreq: 1 uses the band centre only' if n_freq == 1 else f'even with pint_nfreq {n_freq}'
        print(f"WARNING: injected phase of {pulsar_pars['ID']} deviates from PINT by up to {largest*1e9:.0f} ns "
              f"across the band ({hint}); not timing-grade at that level.", flush=True)
    return polycos_path


def write_parfile(pulsar_model, par_path, new_path):
    lines = read_par(par_path)
    par_ephem = par_ephemeris(lines)
    ephem_name, _ = pint_ephemeris(str(pulsar_model.obs.ephem), par_ephem)
    if ephem_name not in (par_ephem, 'builtin'):
        set_par_ephemeris(lines, ephem_name)

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
        # terms in 1/f^2 (Doppler on the DM delay, DM changes, solar wind) and FD-like terms in ln f are fitted
        # exactly; the spline only interpolates what is left, so wide bands need no extra frequencies
        remainder = delays - self.dispersion(nodes)[:, None]
        basis = self.chromatic_basis(nodes, len(nodes))
        self.basis_coef = np.linalg.lstsq(basis, remainder, rcond=None)[0]
        return CubicSpline(np.log(nodes), remainder - basis @ self.basis_coef, axis=0)

    def chromatic_basis(self, freq, n_nodes=None):
        freq = np.asarray(freq, dtype=np.float64)
        n_basis = min(4, n_nodes if n_nodes is not None else self.basis_coef.shape[0])
        log_f = np.log(freq / self.ref_freq)
        columns = [(self.ref_freq / freq)**2 - 1, np.ones_like(log_f), log_f, log_f**2]
        return np.stack(columns[:n_basis], axis=-1)

    def dispersion(self, freq):
        return self.ephemeris_DM * TEMPO_DM_CONST * (1/self.ref_freq**2 - 1/freq**2)

    def delay_table(self, freq):
        freq = np.asarray(freq, dtype=np.float64)
        chromatic = np.broadcast_to(self.dispersion(freq), (len(self.t_grid),) + freq.shape)
        if self.node_delays is not None:
            fitted = np.moveaxis(self.chromatic_basis(freq) @ self.basis_coef, -1, 0)
            chromatic = chromatic + fitted + np.moveaxis(self.node_delays(np.log(freq)), -1, 0)
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
