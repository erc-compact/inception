import numpy as np
import astropy.units as u
from pathlib import Path
import astropy.constants as const
from decimal import Decimal, getcontext
from astropy.time import Time, TimeDelta

from .io_tools import print_exe

getcontext().prec = 40

MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
TEMPO_SITES = {'mk': 'm', 'ao': '3', 'pk': '7', 'jb': '8', 'gb': '1', 'gm': 'r', 'ef': 'g', 'fast': 'k'}


def cheby_nodes(n):
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def cheby_transform(values, n):
    k = np.arange(n)
    T = np.cos(np.outer(np.arange(n), np.pi * (k + 0.5) / n))
    return 2.0 / n * values @ T.T


def cheby_basis(x, n):
    return np.cos(np.outer(np.arccos(np.clip(x, -1, 1)), np.arange(n)))


class PhasePredictor:
    def __init__(self, pulsar_model, psrname, tol=1e-6, edge_tol=1e-4, t2pred_ncoeff_time=16, t2pred_ncoeff_freq=3,
                 polycos_ncoeff=15, min_segment=10.0, max_polycos=450):
        self.pm = pulsar_model
        self.obs = pulsar_model.obs
        self.psrname = psrname
        self.tol = tol
        self.edge_tol = edge_tol
        self.nt = t2pred_ncoeff_time
        self.nf = t2pred_ncoeff_freq
        self.ncoeff = polycos_ncoeff
        self.min_segment = min_segment
        self.max_polycos = max_polycos

        self.prepare_model()
        self.set_ranges()

    def prepare_model(self):
        pm = self.pm
        if not hasattr(pm.binary, 'orbital_delay'):
            pm.binary.orbital_delay = pm.binary.generate_interp()
        self.mode = pm.pulsar_pars.get('mode') or 'python'
        if self.mode == 'pint' and not hasattr(pm, 'polycos'):
            pm.polycos_path = pm.pulsar_pars['polycos']
            pm.get_polyco_interp()

    def set_ranges(self):
        margin = 0.0 if self.mode == 'pint' else 60.0
        self.t_start = -margin
        self.t_end = self.obs.obs_len + margin
        self.f_lo = float(self.obs.low_f)
        self.f_hi = float(self.obs.high_f)
        self.f_ref = float(self.obs.f0)

        first_segment = 1800.0
        if self.pm.binary.period:
            first_segment = min(first_segment, self.pm.binary.period / 4)
        self.first_segment = min(first_segment, self.t_end - self.t_start)

    def DM_delay(self, freq):
        prop = self.pm.prop_effect
        ref = 1/self.obs.high_f**2 if self.pm.pulsar_pars['DM_ref'] == 'top' else 0.0
        return -prop.DM * prop.DM_const * (1/np.asarray(freq, dtype=np.float64)**2 - ref)

    def base_time(self, t_sec):
        if self.pm.pulsar_pars['frame'] == 'bary':
            times = self.obs.obs_start_time_TIME + TimeDelta(t_sec.ravel(), format='sec')
            return self.obs.topo2bary(times, return_mjd=False).reshape(t_sec.shape)
        return t_sec

    def phase(self, t_sec, base, freq):
        delay = self.DM_delay(freq)
        phase = self.pm.get_phase(base + delay)
        if self.mode == 'pint':
            topo_mjd = self.obs.sec2mjd(t_sec.ravel()).reshape(t_sec.shape)
            phase = phase + self.pm.polycos(topo_mjd + delay * u.s.to(u.day))
        return phase

    def mjd_string(self, sec, digits=20):
        return f'{Decimal(self.obs.obs_start) + Decimal(float(sec)) / Decimal(86400):.{digits}f}'

    def fit_t2pred(self, segment):
        n_seg = max(1, int(np.ceil((self.t_end - self.t_start) / segment - 1e-9)))
        edges = self.t_start + (self.t_end - self.t_start) * np.arange(n_seg + 1) / n_seg
        s0, s1 = edges[:-1, None], edges[1:, None]

        xs, ys = cheby_nodes(self.nt), cheby_nodes(self.nf)
        freqs = self.f_lo + (ys + 1) / 2 * (self.f_hi - self.f_lo)
        x_test = np.linspace(-1, 1, 3 * self.nt + 1)
        y_test = np.array([-1.0, 0.0, 1.0])
        f_test = self.f_lo + (y_test + 1) / 2 * (self.f_hi - self.f_lo)

        t_nodes = s0 + (xs[None, :] + 1) / 2 * (s1 - s0)
        t_test = s0 + (x_test[None, :] + 1) / 2 * (s1 - s0)
        t_mid = 0.5 * (s0 + s1)
        times = np.concatenate([t_nodes, t_test, t_mid], axis=1)
        base = self.base_time(times)

        node_ph = np.stack([self.phase(times, base, f)[:, :self.nt] for f in freqs], axis=1)
        test_ph = np.stack([self.phase(times, base, f)[:, self.nt:-1] for f in f_test], axis=1)
        disp_lo = self.phase(times, base, self.f_lo)[:, -1]
        disp_hi = self.phase(times, base, self.f_hi)[:, -1]
        disp = (disp_lo - disp_hi) / (1/self.f_lo**2 - 1/self.f_hi**2)

        values = node_ph - disp[:, None, None] / freqs[None, :, None]**2
        ref = np.floor(values[:, self.nf // 2, self.nt // 2])
        values = values - ref[:, None, None]

        coeffs = np.stack([cheby_transform(cheby_transform(v, self.nt).T, self.nf).T for v in values])
        half = coeffs.copy()
        half[:, 0, :] *= 0.5
        half[:, :, 0] *= 0.5
        model = np.einsum('sji,ti,fj->sft', half, cheby_basis(x_test, self.nt), cheby_basis(y_test, self.nf))
        truth = test_ph - disp[:, None, None] / f_test[None, :, None]**2 - ref[:, None, None]
        err = np.abs(model - truth)

        coeffs[:, 0, 0] += 4 * ref
        return edges, coeffs, disp, float(np.max(err[:, 1])), float(np.max(err[:, [0, 2]]))

    def t2pred_score(self, fit):
        return max(fit[3] / self.tol, fit[4] / self.edge_tol)

    def write_t2pred(self, path):
        segment = self.first_segment
        fit = self.fit_t2pred(segment)
        while self.t2pred_score(fit) > 1 and segment / 2 >= self.min_segment:
            trial = self.fit_t2pred(segment / 2)
            if self.t2pred_score(trial) > 0.75 * self.t2pred_score(fit):
                break
            fit, segment = trial, segment / 2
        edges, coeffs, disp, err_ref, err_edge = fit

        blocks = []
        for k in range(len(coeffs)):
            lines = ['ChebyModel BEGIN', f'PSRNAME {self.psrname}', f'SITENAME {self.obs.telescope_ID.lower()}',
                     f'TIME_RANGE {self.mjd_string(edges[k])} {self.mjd_string(edges[k+1])}',
                     f'FREQ_RANGE {self.f_lo!r} {self.f_hi!r}', f'DISPERSION_CONSTANT {float(disp[k])!r}',
                     f'NCOEFF_TIME {self.nt}', f'NCOEFF_FREQ {self.nf}']
            for ix in range(self.nt):
                lines.append('COEFFS ' + ' '.join(repr(float(coeffs[k, iy, ix])) for iy in range(self.nf)))
            lines.append('ChebyModel END')
            blocks.append('\n'.join(lines))

        with open(path, 'w') as f:
            f.write(f'ChebyModelSet {len(blocks)} segments\n' + '\n'.join(blocks) + '\n')

        print_exe(f'{self.pm.ID} t2pred: {len(blocks)} segments of {segment:.1f} s, max phase error {err_ref:.1e} turns at band centre, {err_edge:.1e} turns at band edges.')
        return path

    def fit_polycos(self, span_min):
        span = span_min * 60.0
        step = 0.9 * span
        n_seg = max(1, int(np.ceil((self.t_end - self.t_start - span) / step - 1e-9)) + 1)

        start = Decimal(self.obs.obs_start)
        tmids = [(start + Decimal(self.t_start + span/2 + k*step) / Decimal(86400)).quantize(Decimal('1e-11')) for k in range(n_seg)]
        mids = np.array([float((tmid - start) * 86400) for tmid in tmids])

        nodes = cheby_nodes(2 * self.ncoeff)
        u_test = np.linspace(-1, 1, 4 * self.ncoeff + 1)
        offsets = np.concatenate([[0.0, -0.5, 0.5], nodes * span / 2, u_test * span / 2])
        times = mids[:, None] + offsets[None, :]
        ph = self.phase(times, self.base_time(times), self.f_ref)

        ph_mid = ph[:, 0]
        f0 = np.round(ph[:, 2] - ph[:, 1], 12)
        n_node = len(nodes)
        dt_node, dt_test = offsets[3:3+n_node], offsets[3+n_node:]
        resid = ph[:, 3:3+n_node] - ph_mid[:, None] - f0[:, None] * dt_node[None, :]
        resid_test = ph[:, 3+n_node:] - ph_mid[:, None] - f0[:, None] * dt_test[None, :]

        fits = np.polynomial.polynomial.polyfit(nodes, resid.T, self.ncoeff - 1).T
        model = np.polynomial.polynomial.polyval(u_test, fits.T)
        err = np.abs(model - resid_test)
        coeffs = fits / (span_min / 2) ** np.arange(self.ncoeff)[None, :]
        rms = np.sqrt(np.mean(err**2, axis=1))
        return tmids, ph_mid, f0, coeffs, rms, float(np.max(err))

    def write_polycos(self, path):
        span_min = int(max(1, min(60, np.floor(self.first_segment / 60))))
        fit = self.fit_polycos(span_min)
        while fit[-1] > self.tol and span_min > 1:
            next_span = max(1, span_min // 2)
            if len(self.fit_polycos_count(next_span)) > self.max_polycos:
                break
            trial = self.fit_polycos(next_span)
            if trial[-1] > 0.75 * fit[-1]:
                break
            fit, span_min = trial, next_span
        tmids, ph_mid, f0, coeffs, rms, err = fit

        turn_offset = max(0, int(-np.floor(np.min(ph_mid))) + 1)
        out = []
        for k, tmid in enumerate(tmids):
            ph_int = int(np.floor(ph_mid[k]))
            frac = ph_mid[k] - ph_int
            frac_r = round(frac, 6)
            if frac_r >= 1.0:
                ph_int, frac, frac_r = ph_int + 1, frac - 1, frac_r - 1
            c = coeffs[k].copy()
            c[0] += frac - frac_r

            t = Time(float(tmid), format='mjd', scale='utc')
            date, hms = t.iso.split()
            yy, mm, dd = date.split('-')
            v_r = self.obs.earth_radial_velocity(float(tmid))[0] if self.pm.pulsar_pars['frame'] == 'bary' else 0.0

            line1 = '{:10.10s} {:>9.9s}{:11.2f}{:20.11f}{:21.6f} {:6.3f}{:7.3f}\n'.format(
                self.psrname, f'{dd}-{MONTHS[int(mm)-1]}-{yy[-2:]}', float(hms.replace(':', '')), float(tmid),
                self.pm.prop_effect.DM, -v_r / const.c.value * 1e4, np.log10(max(rms[k], 1e-12)))
            rphase = '{:13d}'.format(ph_int + turn_offset) + '{:.6f}'.format(frac_r)[1:]
            line2 = '{:20s} {:17.12f}{:>5s}{:5d}{:5d}{:10.3f}{:16s}\n'.format(
                rphase, f0[k], polyco_site(self.obs.tempo_id), span_min, self.ncoeff, self.f_ref, '')
            block = ''
            for i, value in enumerate(c):
                value = 0.0 if abs(value) < 1e-99 else value
                block += '{:25.17e}'.format(value).replace('e', 'D')
                if (i + 1) % 3 == 0:
                    block += '\n'
            if not block.endswith('\n'):
                block += '\n'
            out.append(line1 + line2 + block)

        with open(path, 'w') as f:
            f.write(''.join(out))

        print_exe(f'{self.pm.ID} polycos: {len(tmids)} entries of {span_min} min, max phase error {err:.1e} turns at {self.f_ref:.3f} MHz.')
        return path

    def fit_polycos_count(self, span_min):
        span = span_min * 60.0
        n_seg = max(1, int(np.ceil((self.t_end - self.t_start - span) / (0.9 * span) - 1e-9)) + 1)
        return range(n_seg)


def polyco_site(tempo_id):
    code = TEMPO_SITES.get(tempo_id, tempo_id)
    return str(ord(code) - ord('a') + 10) if len(code) == 1 and code.isalpha() else code


def predictor_name(pulsar_model, parfile_path):
    if parfile_path and Path(parfile_path).suffix == '.par' and Path(parfile_path).exists():
        with open(parfile_path) as f:
            for line in f:
                columns = line.split()
                if len(columns) > 1 and columns[0] in ('PSR', 'PSRJ', 'PSRB'):
                    return columns[1].lstrip('J')
    return str(pulsar_model.ID)


def create_predictor(pulsar_model, kind, output_path, parfile_path=None):
    predictor = PhasePredictor(pulsar_model, predictor_name(pulsar_model, parfile_path))
    path = f'{output_path}/{pulsar_model.ID}_predictor.{kind}'
    if kind == 't2pred':
        return predictor.write_t2pred(path)
    return predictor.write_polycos(path)
