import os
import sys
import json
import argparse
import subprocess

class FoldExec:
    def __init__(self) -> None:
        pass


    """
    psrfold_fil2 --dmboost 250 --plotx -v -t 12 --candfile {} -n {} {} {} --template {} --clfd 8 -L {} --fillPatch rand -f {} --rfi zdot {} --fd {} --td {}".format(
                    pred_file, nsubband, nbins_string, beam_tag, TEMPLATE, subint_length, input_filenames, zap_string, fscrunch, tscrunch)
    """