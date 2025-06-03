import argparse
from setup_manager import SetupManager
from signal_injector import InjectSignal
from pulsar_par_parser import PulsarParParser


if __name__=='__main__':
    pulsar_params = PulsarParParser({'ID': 'help', 'P0': 1, 'SNR': 1}, {}).parser.print_help()

    parser = argparse.ArgumentParser(prog='Inception injector',
                                     epilog=f'Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de.')
    parser.add_argument('--signal', metavar='file', required=True, help='text file containing pulsar parameters to inject')
    parser.add_argument('--fb', metavar='file', required=True, help='path to filterbank where signal is injected')
    parser.add_argument('--output', metavar='directory', required=True, help='output directory for injected filterbank')
    parser.add_argument('--ephem', metavar='file', required=False, default='builtin', help='JPL ephemeris file for solar system (.bsp)')
    parser.add_argument('--ncpu', metavar='integer', required=False, default=1, type=int, help='number of cpus')
    args = parser.parse_args()

    setup = SetupManager(args.signal, args.fb, args.ephem, args.output)
   
    injector = InjectSignal(setup, args.ncpu)
    injector.parallel_inject()
    injector.combine_files()
    

    
