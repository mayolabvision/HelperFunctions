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

echo "Job started at $(date)"

module purge
module load matlab/R2023a

# Define the varargin parameters
RAW_PATH='/ix1/pmayo/lab_NHPdata'
OUT_PATH='/ix1/pmayo/OneDrive/DATA'
CSV_PATH='/ix1/pmayo/lab_NHPdata/RECORDING_INFO.csv'
NET_PATH='/ihome/pmayo/knoneman/Packages/nasnet'
NEV_PATH='/ihome/pmayo/knoneman/Packages/nevutils'

matlab -nodisplay -r "try, process_RippleRecording_KKN('$1', '$2', '$3', 'RAW_DATA_PATH', '$RAW_PATH', 'OUT_DATA_PATH', '$OUT_PATH', 'RECD_CSV_PATH', '$CSV_PATH', 'NASNET_PATH', '$NET_PATH', 'NEVUTIL_PATH', '$NEV_PATH'); catch, exit(1); end; exit;"

echo "DONE"
