#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=mtlab
#SBATCH --error=/ix1/pmayo/matlab/outfiles/error_%A.err
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-00:09:59

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."

module purge
module load matlab/R2023a

# Define the varargin parameters
RAW_PATH='/ix1/pmayo/lab_NHPdata'
OUT_PATH='/ix1/pmayo/OneDrive/DATA'
CSV_PATH='/ix1/pmayo/lab_NHPdata/RECORDING_INFO.csv'
NET_PATH='/ihome/pmayo/knoneman/Packages/nasnet/networks'

PATH1='/ihome/pmayo/knoneman/Packages/HelperFunctions'
PATH2='/ihome/pmayo/knoneman/Packages'

matlab -nodisplay -r "try, process_RippleRecording_KKN('$1', '$2', '$3', 'RAW_PATH', '$RAW_PATH', 'OUT_PATH', '$OUT_PATH', 'CSV_PATH', '$CSV_PATH', 'NET_PATH', '$NET_PATH', 'PATH1', '$PATH1', 'PATH2', '$PATH2'); catch, exit(1); end; exit;"

echo "DONE"
