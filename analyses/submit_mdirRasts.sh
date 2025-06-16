#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=mdir
#SBATCH --error=/ix1/pmayo/matlab/outfiles/out_%A_%a.out
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-01:59:59
#SBATCH --array=0-99

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "My SLURM_ARRAY_TASK_COUNT is $SLURM_ARRAY_TASK_COUNT"
echo "Job started at $(date)"

module purge
module load matlab/R2023a

#########INPUTS##########
SESSION="$1"
IMEC="${2:-0}"
RUN_TYPE="${3:-unleashed}"
SWEEP_NAME="${4:-none}"

HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

echo "SESSION: $SESSION"
echo "IMEC: $IMEC"
echo "RUN_TYPE: $RUN_TYPE"

if [[ "$RUN_TYPE" == "sweep" ]]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/tables/${SESSION}_${RUN_TYPE}_${SWEEP_NAME}.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/figs/kilosort4_${RUN_TYPE}/${SWEEP_NAME}/mdir_rasters"
else
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/tables/${SESSION}_${RUN_TYPE}.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/figs/kilosort4_${RUN_TYPE}/mdir_rasters"
fi

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"

########################

ALIGN="stim"
echo "ALIGN: $ALIGN"
# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
addpath(genpath('$HELPERS_PATH'));
fprintf('Running ia_mdirRasters for $1\n');
ia_mdirRasters('$DATA_PATH', ...
    'IMEC', $IMEC, ...
    'ALIGN', '$ALIGN', ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

ALIGN="sacc"
echo "ALIGN: $ALIGN"
matlab -nodisplay <<EOF
addpath(genpath('$HELPERS_PATH'));
fprintf('Running ia_mdirRasters for $1\n');
ia_mdirRasters('$DATA_PATH', ...
    'IMEC', $IMEC, ...
    'ALIGN', '$ALIGN', ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

echo "DONE"
