import os
import sys
import json
import argparse
import subprocess
import numpy as np
from pathlib import Path
from collections import namedtuple

from injector.io_tools import merge_filterbanks


class InjectorSetup:
    def __init__(self, search_args, inject_file, xml_dir, work_dir):
        self.out = work_dir
        self.xml_dir = xml_dir

        args = self.parse_JSON(search_args)
        self.inject_file = self.parse_JSON(inject_file)
        self.processing_args, self.ID  = args['processing_args'], args['processing_id']
        self.beam_data = self.get_beam(args['data']['pointings'])
        
        self.injection_ID = self.inject_file['injection_id']
        self.data_ID = args['processing_id']

        self.create_XML_list(args['data']['n_cbeams'], args['data']['n_ibeams'])

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

    def get_beam(self, pointings):
        beams = []
        for pointing in pointings:
            for beam in pointing['beams']:
                beam_data = [pointing['id'], beam['name']]
                for filterbank in beam['data_products']:
                    beam_data.append(filterbank['filename'])
                beams.append(beam_data)

        selected_beam = np.random.choice(beams)
        return selected_beam
    
    def create_DDplan(self):
        DMRange = namedtuple("DMRange", ["low_dm", "high_dm", "dm_step", "tscrunch"])

        segments = []
        plan = self.processing_args['ddplan']
        for line in plan.splitlines():
            low_dm, high_dm, dm_step, tscrunch = list(map(float, line.split()[:4]))
            segments.append(DMRange(low_dm, high_dm, dm_step, tscrunch))

        return list(sorted(segments, key=lambda x: x.tscrunch))

    def create_XML_list(self, n_cbeams, n_ibeams):
        ddplan = self.create_DDplan()
        xml_file_names = [f'overview_dm_{dm_range.low_dm:.6f}_{dm_range.high_dm:.6f}.xml' for dm_range in ddplan]

        pointing_id, inj_beam_name, *_ = self.beam_data
        xml_file_paths = []
        for nbeam in range(n_cbeams):
            for xml_name in xml_file_names:
                beam_i = f'cfbf{nbeam:05g}'
                if beam_i != inj_beam_name:
                    xml_file_paths.append(f'{self.xml_dir}/{pointing_id}/XML_FILES/{beam_i}/{xml_name}')
        if n_ibeams:
            for xml_name in xml_file_names:
                xml_file_paths.append(f'{self.xml_dir}/{pointing_id}/XML_FILES/ifbf00000/{xml_name}')

        for xml_name in xml_file_names:
            xml_file_paths.append(f'{self.out}/{self.injection_ID}_{self.data_ID}_{xml_name}')

        outfile = f'{self.out}/candidates_{self.injection_ID}_{self.data_ID}.ascii'
        np.savetxt(outfile, xml_file_paths)

    def merge_data_products(self):
        pointing_id, inj_beam_name, *fb_names = self.beam_data
        output_name = f'{self.out}/{pointing_id}_{inj_beam_name}_{self.injection_ID}_{self.data_ID}.fil'
        ### change fb_names path with output dir Path(...)
        merge_filterbanks(fb_names, output_name)
        return output_name

    def run_injector(self, ephem, ncpus):
        fb_path = self.merge_data_products()

        script_path = './inception/injector'
        inputs = f"--signal={self.inject_file} --fb={fb_path} --ephem={ephem} --output={self.out} --ncpu={ncpus}"
        cmd = f"python3 {script_path}/SCRIPT_inject_pulsars.py {inputs}"

        subprocess.run(cmd, shell=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='offline injection pipeline setup',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--search_args', metavar='file', required=True, help='JSON file with search parameters')
    parser.add_argument('--injection_file', metavar='file', required=True, help='JSON file with injection plan')
    parser.add_argument('--xml_dir', metavar='dir', required=True, help='directory of xml files')
    parser.add_argument('--work_dir', metavar='dir', required=True, help='work directory')
    parser.add_argument('--ephem', metavar='file', required=False, default='builtin', help='solar system ephemeris file')
    parser.add_argument('--ncpus', metavar='int', required=False, default=1, type=int, help='number of threads for injection')
    args = parser.parse_args()

    inj_setup = InjectorSetup(args.search_args, args.injection_file, args.data_dir, args.work_dir)
    inj_setup.run_injector(args.ephem, args.ncpus)

