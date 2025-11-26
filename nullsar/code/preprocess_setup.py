import sys, os
import argparse
from collections import defaultdict

from nullsar_tools import parse_JSON


class Setup:
    def __init__(self, processing_args, filterbanks, out_dir, work_dir):

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir

        self.processing_args_path = processing_args
        self.processing_args = parse_JSON(processing_args)

        self.fb = self.parse_fb(filterbanks)

    @staticmethod
    def parse_fb(filterbanks):
        parsed = defaultdict(list)

        with open(filterbanks) as f:
            for line in f:
                tag, path = line.strip().split(maxsplit=1)
                parsed[tag].append(path)

        return parsed

    def create_dir(self):
        os.makedirs(f'{self.out_dir}/PROCESSING', exist_ok=True)
        os.makedirs(f'{self.out_dir}/SETUP', exist_ok=True)

        for tag in self.fb.keys():
            os.makedirs(f'{self.out_dir}/PROCESSING/{tag}', exist_ok=True)
            os.makedirs(f'{self.out_dir}/PROCESSING/{tag}/01_FILES', exist_ok=True)

    def create_tag_file(self):
        for loc in [self.work_dir, f'{self.out_dir}/SETUP']:
            with open(f'{loc}/tags.txt', 'w') as f:
                for tag in sorted(self.fb.keys()):
                    f.write(tag + "\n")

    def create_fb_list(self):
        for tag in self.fb.keys():
            with open(f'{self.out_dir}/PROCESSING/{tag}/01_FILES/files.txt', 'w') as f:
                for fb in self.fb[tag]:
                    f.write(fb + "\n")

    def run_setup(self):

        self.create_dir()

        self.create_tag_file()

        self.create_fb_list()



if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='NULLSAR - preprocessing setup',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')
    parser.add_argument('--processing_args', metavar='file', required=True, help='JSON file with preprocessing parameters')
    parser.add_argument('--filterbanks', metavar='file', required=True,  help='.txt file containg filterbanks and processing tags')
    parser.add_argument('--out_dir', metavar='dir', required=False, default='cwd', help='output directory')
    parser.add_argument('--work_dir', metavar='dir', required=False, default='cwd', help='work directory')

    args = parser.parse_args()
    setup_exec = Setup(args.processing_args, args.filterbanks, args.out_dir, args.work_dir)

    setup_exec.run_setup()