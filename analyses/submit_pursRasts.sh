#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=purs
#SBATCH --error=/ix1/pmayo/matlab/outfiles/out_%A_%a.out
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
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
IMEC="${2:-0}"
RUN_TYPE="${3:-unleashed}"
SWEEP_NAME="${4:-none}"

ALIGN="targ"
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

echo "SESSION: $SESSION"
echo "IMEC: $IMEC"
echo "RUN_TYPE: $RUN_TYPE"

if [[ "$RUN_TYPE" == "sweep" ]]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/tables/${SESSION}_${RUN_TYPE}_${SWEEP_NAME}.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/figs/kilosort4_${RUN_TYPE}/${SWEEP_NAME}/purs_rasters"
else
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/tables/${SESSION}_${RUN_TYPE}.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/figs/kilosort4_${RUN_TYPE}/purs_rasters"
fi

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"

########################

PURE_ONLY="0"
echo "ALL TRIALS"
matlab -nodisplay <<EOF
addpath(genpath('$HELPERS_PATH'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'IMEC', $IMEC, ...
    'ALIGN', '$ALIGN', ...
    'PURE_ONLY', logical(str2double('$PURE_ONLY')), ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

PURE_ONLY="1"
echo "PURE ONLY"
matlab -nodisplay <<EOF
addpath(genpath('$HELPERS_PATH'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'IMEC', $IMEC, ...
    'ALIGN', '$ALIGN', ...
    'PURE_ONLY', logical(str2double('$PURE_ONLY')), ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')), ...
    'N_CHUNKS', str2double(getenv('SLURM_ARRAY_TASK_COUNT')));
exit
EOF

echo "DONE"
