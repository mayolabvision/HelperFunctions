#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --job-name=rfmp
#SBATCH --error=/ix1/pmayo/matlab/outfiles/error_%A_%a.err
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-00:59:59
#SBATCH --array=0-49

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "My SLURM_ARRAY_TASK_COUNT is $SLURM_ARRAY_TASK_COUNT"
echo "Job started at $(date)"

module purge
module load matlab/R2023a

#########INPUTS##########

SESSION="$1"
IMEC="${2:-0}"
CORRECTED="${3:-0}"

echo "SESSION: $SESSION"
echo "IMEC: $IMEC"
echo "CORRECTED: $CORRECTED"

########################

# Define the varargin parameters
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

if [ "$CORRECTED" -eq 0 ]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_nodrift.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort_nodrift/rfmp/unit_heatmaps"
elif [ "$CORRECTED" -eq 1 ]; then
    DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}_medicine.mat"
    FIG_PATH="/ix1/pmayo/lab_NHPdata/${1}/figs/kilosort_medicine/rfmp/unit_heatmaps"
else
    echo "Error: Invalid CORRECTED value '$CORRECTED'. Must be 0 or 1."
    exit 1
fi

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"

matlab -nodisplay <<EOF
try
    addpath(genpath('$HELPERS_PATH'));
    fprintf('Running ia_rfMaps for $1\n');
    ia_rfMaps('$DATA_PATH', ...
        'IMEC', $IMEC, ...
        'FIG_PATH', '$FIG_PATH', ...
        'JOB_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')));
catch err
    disp('ERROR in ia_rfMaps:');
    disp(getReport(err));
    exit(1);
end
exit
EOF

echo "DONE"
