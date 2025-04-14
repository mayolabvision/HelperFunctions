#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --job-name=npy2mat
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --error=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-01:59:59

echo "Job started at $(date)"

module purge
conda activate /ihome/pmayo/knoneman/.conda/envs/npy2mat

# Define the varargin parameters
KILO4_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_imec${2}/kilosort4"

# Convert .npy files in kilosort4 directory to .mat
python -c "import sys; sys.path.append('/ihome/pmayo/knoneman/Packages/HelperFunctions/utils'); from convert_npy_to_mat import convert_npy_to_mat; convert_npy_to_mat('$KILO4_PATH')"

echo "DONE"
