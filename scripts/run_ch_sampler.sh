#!/bin/bash                                                                                                                            
#
#SBATCH --job-name=master_mc_ch
#SBATCH --output=master_mc_ch.out
#SBATCH --mail-type=END,FAIL
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=8
# SBATCH --gres=gpu:1
# SBATCH -C k40
#SBATCH --time=1:00:00
#SBATCH --mem=8gb

output_dir=${HOME}/maDGiCart-CH/scripts/output3D

for i in {1..128}
do
    work_dir=${output_dir}/run_${i}
    mkdir ${work_dir}
    cp base_ch_script.sh ${work_dir}/run_${i}.sh
    sbatch ${work_dir}/run_${i}.sh ${work_dir}
done
