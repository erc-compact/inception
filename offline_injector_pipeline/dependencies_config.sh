#!/bin/bash
tmp_dir=$1

python3.6 -m pip install --target $tmp_dir sympy
python3.6 -m pip install --target $tmp_dir astropy
python3.6 -m pip install --target $tmp_dir pandas
export PYTHONPATH=$tmp_dir:$PYTHONPATH