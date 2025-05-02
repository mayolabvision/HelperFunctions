#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=rfmp
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
IMEC="${2:-0}"
RUN_TYPE="${3:-unleashed}"

HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

echo "SESSION: $SESSION"
echo "IMEC: $IMEC"
echo "RUN_TYPE: $RUN_TYPE"

DATA_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/${SESSION}_${RUN_TYPE}.mat"
FIG_PATH="/ix1/pmayo/lab_NHPdata/${SESSION}/figs/kilosort4_${RUN_TYPE}/rfmp_heatmaps"

echo "DATA_PATH: $DATA_PATH"
echo "FIG_PATH: $FIG_PATH"

########################

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
