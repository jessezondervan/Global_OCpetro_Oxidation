#!/bin/bash 
#SBATCH --job-name=Global_OCpetro_oxidation_model
#SBATCH --partition=short
#SBATCH --time=00:45:00 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48

cd $TMPDIR || exit 1

rsync -av $DATA/input_global ./

module purge
module load Anaconda3//2021.11

source activate $DATA/carbon-env-newer

python3 input_global/Glob_newmethod_parr_globalresidual.py

rsync -av --exclude=input_global --exclude='glob_petrodenud.tif' --exclude='glob_petrodenud_wo_alluv.tif' ./ $DATA
