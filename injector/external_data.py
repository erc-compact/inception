import os
import sys
import warnings
import astropy.units as u
from astropy.time import Time
from astropy.utils import iers
from astropy.coordinates import EarthLocation
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.data import conf as data_conf, CacheMissingWarning, is_url_in_cache, download_file

from .io_tools import print_exe


TEMPO2_SITES = {'Meerkat': (5109360.133, 2006852.586, -3238948.127),
                'Arecibo': (2390487.080, -5564731.357, 1994720.633),
                'Parkes': (-4554231.5, 2816759.1, -3454036.3),
                'Jodrell': (3822626.04, -154105.65, 5086486.04),
                'GBT': (882589.289, -4924872.368, 3943729.418),
                'GMRT': (1656342.30, 5797947.77, 2073243.16),
                'Effelsberg': (4033949.5, 486989.4, 4900430.8),
                'FAST': (-1668557.2, 5506838.5, 2744934.9)}

_state = None


def configure(cache=None, offline=False):
    if cache:
        os.environ['INJECTOR_CACHE'] = os.path.abspath(cache)
    if offline:
        os.environ['INJECTOR_OFFLINE'] = '1'
    apply()


def apply():
    global _state
    cache = os.environ.get('INJECTOR_CACHE', '')
    offline = os.environ.get('INJECTOR_OFFLINE') == '1'
    if _state == (cache, offline):
        return

    if cache:
        os.makedirs(os.path.join(cache, 'astropy'), exist_ok=True)
        os.environ['XDG_CACHE_HOME'] = cache
    else:
        warnings.filterwarnings('ignore', category=CacheMissingWarning)

    if not offline and not online_tables():
        print_exe('astropy downloads failed, using local IERS tables')
        os.environ['INJECTOR_OFFLINE'] = '1'
        offline = True
    if offline:
        offline_tables()

    _state = (cache, offline)


def online_tables():
    data_conf.remote_timeout = 5
    iers.conf.iers_degraded_accuracy = 'warn'
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', AstropyWarning)
            warnings.simplefilter('error', iers.IERSWarning)
            iers.IERS_Auto.open().ut1_utc(Time.now())
    except Exception:
        return False
    return True


def offline_tables():
    data_conf.allow_internet = False
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    iers.conf.iers_degraded_accuracy = 'ignore'
    iers.earth_orientation_table.set(local_table())


def local_table():
    for url in (iers.conf.iers_auto_url, iers.conf.iers_auto_url_mirror):
        if is_url_in_cache(url):
            try:
                return iers.IERS_A.open(download_file(url, cache=True))
            except Exception:
                pass
    return iers.IERS_B.open()


def get_observatory(name):
    if name in TEMPO2_SITES:
        return EarthLocation.from_geocentric(*TEMPO2_SITES[name], unit=u.m)
    try:
        return EarthLocation.of_site(name)
    except Exception as err:
        sys.exit(f'Telescope {name} not in the built in table or the astropy site registry: {err}')
