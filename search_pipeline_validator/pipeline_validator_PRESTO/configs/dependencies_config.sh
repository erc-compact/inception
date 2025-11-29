#!/bin/bash
tmp_dir=$1
python_cmd=$2

$python_cmd -m pip install --timeout 100 --retries 5 --target "$tmp_dir" sympy
$python_cmd -m pip install --timeout 100 --retries 5 --target "$tmp_dir" astropy
$python_cmd -m pip install --timeout 100 --retries 5 --target "$tmp_dir" pandas
export PYTHONPATH="$tmp_dir:${PYTHONPATH:-}"


