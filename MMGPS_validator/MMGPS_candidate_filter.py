import os
import re
import argparse
import subprocess
import pandas as pd
from pathlib import Path

from ..search_pipeline_validator import pipeline_tools as inj_tools


class CandidateFilterProcess:
    def __init__(self,  processing_args, out_dir, work_dir, injection_number):
        self.processing_args_path = processing_args
        self.processing_args = inj_tools.parse_JSON(processing_args)

        self.out_dir = os.getcwd() if out_dir == 'cwd' else out_dir
        self.work_dir = os.getcwd() if work_dir == 'cwd' else work_dir
        
        self.injection_number = injection_number