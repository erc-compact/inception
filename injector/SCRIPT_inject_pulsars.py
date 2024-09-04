import sys
from pathlib import Path
sys.path.insert(0, Path(__file__).parent.parent)

import getopt
from injector.setup_manager import SetupManager
from injector.signal_injector import InjectSignal


if __name__=='__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["signal=", "filterbank=", "ephem=", "output=", 'ncpu='])
    except Exception as err:
        sys.exit(err)

    ad = dict(opts) 
    if not ad.get('-h', 1):
        sys.exit('If you need help come to my office (E0.04). :D') 
        
    setup = SetupManager(ad['--signal'], ad['--filterbank'],  ad.get('--ephem', 'builtin'),  ad['--output'])
   
    injector = InjectSignal(setup, ad.get('--ncpu', 1))
    injector.parallel_inject()
    injector.combine_files()
    

    
