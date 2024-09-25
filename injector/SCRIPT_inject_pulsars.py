import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

import argparse
from injector.setup_manager import SetupManager
from injector.signal_injector import InjectSignal


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Inception injector',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--signal', metavar='file', required=True, help='text file containing pulsar parameters to inject')
    parser.add_argument('--filterbank', metavar='file', required=True, help='path to filterbank where signal is injected')
    parser.add_argument('--output', metavar='directory', required=True, help='output directory for injected filterbank')
    parser.add_argument('--ephem', metavar='file', required=False, default='builtin', help='JPL ephemeris file for solar system (.bsp)')
    parser.add_argument('--ncpu', metavar='integer', required=False, default=1, type=int, help='number of cpus')
    args = parser.parse_args()

    setup = SetupManager(args.signal, args.filterbank, args.ephem, args.output)
   
    injector = InjectSignal(setup, args.ncpu)
    injector.parallel_inject()
    injector.combine_files()
    

    
