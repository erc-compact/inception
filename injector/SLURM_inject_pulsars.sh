#!/bin/bash
#SBATCH --partition=short.q
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=7G
#SBATCH --time=4:00:00
#SBATCH --job-name=inject

fb_path="/path/filterbank.fil"
inj_path="/path/example.inject"
ephem_path="/path/de440.bsp"
profile_path="/path/profile.npy"

sing_img="/path/singularity.sif"

tmp_dir="/tmp"

rm -rf $tmp_dir
mkdir -p $tmp_dir

rsync -Pav $fb_path $tmp_dir
rsync -Pav $inj_path $tmp_dir
rsync -Pav $ephem_path $tmp_dir
rsync -Pav $profile_path $tmp_dir


command="python3 SCRIPT_inject_pulsars.py"  
inputs="--signal=$tmp_dir/example.inject --filterbank=$tmp_dir/filterbank.fil --ephem=$tmp_dir/de440.bsp --output=$tmp_dir --ncpu=$SLURM_NTASKS"
singularity exec -H $HOME:/home -B $sing_img $command $inputs

# warning: I forgot to delete temporary files created during the injection. Ignore .tmpfil files
output_dir="/path/output"
rsync -Pav  $tmp_dir/*.par  $output_dir
rsync -Pav  $tmp_dir/filterbank_pulsar1.fil  $output_dir

rm -rf $tmp_dir
