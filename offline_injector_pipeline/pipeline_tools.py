import sys
import json
from collections import namedtuple

class PipelineTools:
    def __init__(self, search_args):
        self.parse_search_args(search_args)
        
    def parse_search_args(self, search_args):
        self.search_args = search_args
        args = self.parse_JSON(search_args)
        self.data_ID = args['processing_id']
        self.processing_args = args['processing_args']
        self.fold_args = args['fold_args']
        self.multi_beam = args['multi_beam_args']
        self.data = args['data']

    def parse_JSON(self, json_file):
        try:
            with open(json_file, 'r') as file:
                pars = json.load(file)
        except FileNotFoundError:
            sys.exit(f'Unable to find {json_file}.')
        except json.JSONDecodeError:
            sys.exit(f'Unable to parse {json_file} using JSON.')
        else:
            return pars
        
    def create_DDplan(self):
        DMRange = namedtuple("DMRange", ["low_dm", "high_dm", "dm_step", "tscrunch"])

        segments = []
        plan = self.processing_args['ddplan']
        for line in plan.splitlines():
            low_dm, high_dm, dm_step, tscrunch = list(map(float, line.split()[:4]))
            segments.append(DMRange(low_dm, high_dm, dm_step, tscrunch))

        return list(sorted(segments, key=lambda x: x.tscrunch))
    
