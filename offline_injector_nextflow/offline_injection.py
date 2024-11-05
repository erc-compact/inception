import os
import sys
import json
import argparse
import subprocess


# choose random beam, 
# get xml_file list
# merge filterbanks
# run injector

class InjectorSetup:
    def __init__(self, search_args):
        args = self.parse_JSON(search_args)
        self.pointings = args['data']['pointings']


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