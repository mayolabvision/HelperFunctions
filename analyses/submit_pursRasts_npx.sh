#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=purs
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
ALIGN="${2:-targ}"
IMEC="${3:-0}"
DRIFT_CORRECT_TYPE="${4:-None}"

echo "SESSION: $SESSION"
echo "IMEC: $IMEC"
echo "ALIGN: $ALIGN"
echo "DRIFT_CORRECT_TYPE: $DRIFT_CORRECT_TYPE"

########################

# Define the varargin parameters
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

if [ "$DRIFT_CORRECT_TYPE" == "medicine" ]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_medicine.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort4_medicine/purs/unit_rasters"
elif [ "$DRIFT_CORRECT_TYPE" == "kilosort" ]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_kilosort.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort4_kilosort/purs/unit_rasters"
else
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_none.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort4_none/purs/unit_rasters"
fi

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"

matlab -nodisplay <<EOF
addpath(genpath('$HELPERS_PATH'));
fprintf('Running ia_pursRasters for $1\n');
ia_pursRasters('$DATA_PATH', ...
    'IMEC', $IMEC, ...
    'ALIGN', '$ALIGN', ...
    'FIG_PATH', '$FIG_PATH', ...
    'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')));
exit
EOF

echo "DONE"
