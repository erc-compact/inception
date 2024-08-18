import sys
import getopt

from injector.setup_manager import SetupManager
from injector.signal_injector import InjectSignal


if __name__=='__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["signal=", "filterbank=", "ephem=", "output=", 'ncpu='])
    except:
        sys.exit(sys.exc_info()[0])

    ad = dict(opts)    
    setup = SetupManager(ad['--signal'], ad['--filterbank'], ad['--ephem'],  ad['--output'])
   
    injector = InjectSignal(setup, ad['--ncpu'])
    injector.parallel_inject()
    injector.combine_files()
    

    