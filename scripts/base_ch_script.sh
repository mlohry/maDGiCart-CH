#!/bin/bash                                                                                                                            
#
#SBATCH --job-name=mc_ch
#SBATCH --output=mc_ch.out
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:1
#SBATCH -C k40
#SBATCH --time=01:00:00
#SBATCH --mem=32gb

driver_dir=${HOME}/maDGiCart-CH/scripts
build_dir=${HOME}/maDGiCart-CH-build/
sif_dir=${HOME}/maDGiCart-CH

singularity exec --nv --env LC_ALL=C ${sif_dir}/myimage.sif python3 -u ${driver_dir}/sampler.py $1 ${build_dir}
