#!/bin/bash
tmp_dir=$1

$2 -m pip install --target $tmp_dir sympy
$2 -m pip install --target $tmp_dir astropy
$2 -m pip install --target $tmp_dir pandas
export PYTHONPATH="$tmp_dir:${PYTHONPATH:-}"
