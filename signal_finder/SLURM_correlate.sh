#!/bin/bash
#SBATCH --partition=short.q
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=7G
#SBATCH --time=4:00:00
#SBATCH --job-name=correlate

fb_path="/path/filterbank.fil"
sing_img="/path/singularity.sif"

tmp_dir="/tmp"

rm -rf $tmp_dir
mkdir -p $tmp_dir

rsync -Pav $fb_path $tmp_dir

command="python3 SCRIPT_correlate.py"  
inputs="--filterbank=$tmp_dir/filterbank.fil --output=$tmp_dir --ncpu=$SLURM_NTASKS --DM_min=0 --DM_max=500 --DM_res=2"
singularity exec -H $HOME:/home -B $sing_img $command $inputs


output_dir="/path/output"
rsync -Pav  $tmp_dir/stacked_DM.png  $output_dir
rsync -Pav  $tmp_dir/stacked_DM.npy  $output_dir

rm -rf $tmp_dir
