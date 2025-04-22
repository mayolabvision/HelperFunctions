#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=mdir
#SBATCH --error=/ix1/pmayo/matlab/outfiles/error_%A_%a.err
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-00:59:59
#SBATCH --array=0-99

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "My SLURM_ARRAY_TASK_COUNT is $SLURM_ARRAY_TASK_COUNT"
echo "Job started at $(date)"

module purge
module load matlab/R2023a

#########INPUTS##########

SESSION="$1"
ALIGN="$2"
IMEC="${3:-0}"
CORRECTED="${4:-0}"

echo "SESSION: $SESSION"
echo "IMEC: $IMEC"
echo "ALIGN: $ALIGN"
echo "CORRECTED: $CORRECTED"

########################

# Define the varargin parameters
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

if [ "$CORRECTED" -eq 0 ]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_baseline.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort_baseline/mdir/unit_rasters"
elif [ "$CORRECTED" -eq 1 ]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_medicine.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort_medicine/rfmp/unit_rasters"
else
    echo "Error: Invalid CORRECTED value '$CORRECTED'. Must be 0 or 1."
    exit 1
fi

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"

# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
addpath(genpath('$HELPERS_PATH'));
fprintf('Running ia_mdirRasters for $1\n');
ia_mdirRasters('$DATA_PATH', ...
    'IMEC', $IMEC, ...
    'ALIGN', '$ALIGN', ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')));
exit
EOF

echo "DONE"
