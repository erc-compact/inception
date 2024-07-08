#!/bin/bash
#SBATCH --partition=short.q
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=10G
#SBATCH --time=4:00:00
#SBATCH --job-name=psr

fb_name="controlfb.fil"
fb_path="/hercules/scratch/rsenzel/signal_inject/$fb_name"
par_path="/u/rsenzel/pulsar_inject/test1.pulsar"

sing_img="/hercules/:/hercules/ /hercules/scratch/vishnu/singularity_images/pulsar-miner_turing-sm75.sif"
tmp_in="/tmp/rsenzel/input_data"
tmp_out="/tmp/rsenzel/output_data" 

rm -rf $tmp_in
rm -rf $tmp_out
mkdir -p $tmp_in
mkdir -p $tmp_out
rsync -Pav $fb_path $tmp_in
rsync -Pav $par_path $tmp_in

command="python3 script_inject_signal.py" 
inputs="--signal $tmp_in/example.pulsar --filterbank=$tmp_in/$fb_name --output=$tmp_out"
singularity exec -H $HOME:/home -B $sing_img $command $inputs

output_dir="/hercules/scratch/rsenzel/signal_inject"
mkdir -p $output_dir
rsync -Pav  $tmp_out/*.fil  $output_dir
rm -rf $tmp_in
rm -rf $tmp_out
