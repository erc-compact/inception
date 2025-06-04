import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

from injector.setup_manager import SetupManager
from injector.signal_injector import InjectSignal
from injector.pulsar_par_parser import PulsarParParser


if __name__=='__main__':
    pulsar_params = PulsarParParser({'ID': 'help', 'P0': 1, 'SNR': 1}, {}).parser.print_help()

    parser = argparse.ArgumentParser(prog='Inception injector',
                                     epilog=f'Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de.')
    parser.add_argument('--signal', metavar='file', required=True, help='text file containing pulsar parameters to inject')
    parser.add_argument('--fb', metavar='file', required=True, help='path to filterbank where signal is injected')
    parser.add_argument('--output', metavar='directory', required=True, help='output directory for injected filterbank')
    parser.add_argument('--ephem', metavar='file', required=False, default='builtin', help='JPL ephemeris file for solar system (.bsp)')
    parser.add_argument('--ncpu', metavar='integer', required=False, default=1, type=int, help='number of cpus')

    parser.add_argument('--gulp_size_GB', metavar='float', required=False, default=0.1, type=int, help='injection gulp size in GB')
    parser.add_argument('--stats_samples', metavar='integer', required=False, default=1e6, type=float, help='number of samples to use for fb statistics')
    args = parser.parse_args()

    setup = SetupManager(args.signal, args.fb, args.ephem, args.output)
   
    injector = InjectSignal(setup, args.ncpu, args.gulp_size_GB, args.stats_samples)
    injector.parallel_inject()
    injector.combine_files()
    

    
