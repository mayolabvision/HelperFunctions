#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=mtlab
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --error=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-02:00:00

echo "Job started at $(date)"

module purge
module load matlab/R2023a
conda activate /ihome/pmayo/knoneman/.conda/envs/npy2mat

# Assign variables from positional parameters
SESSION="$1"
RUN_TYPE="${2:-behavOnly}"

echo "SESSION: $SESSION"
echo "RUN_TYPE: $RUN_TYPE"

# Define the varargin parameters
RAW_PATH='/ix1/pmayo/lab_NHPdata'
OUT_PATH='/ix1/pmayo/lab_NHPdata'
NEV_PATH='/ihome/pmayo/knoneman/Packages/nevutils'
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'

##########################################################################3

# First MATLAB call: run process_fullRecording
matlab -nodisplay <<EOF
try
    addpath(genpath('$HELPERS_PATH'));
    fprintf('Running process_fullRecording for $1\n');
    process_fullRecording('${SESSION}', ...
        'RAW_DATA_PATH', '$RAW_PATH', ...
        'OUT_DATA_PATH', '$OUT_PATH', ...
        'NEVUTIL_PATH', '$NEV_PATH', ...
        'RUN_TYPE', '$RUN_TYPE');
catch err
    disp('ERROR in process_fullRecording:');
    disp(getReport(err));
    exit(1);
end
exit
EOF

matlab -nodisplay -r "try, addpath(genpath('$HELPERS_PATH')); ia_trialOutcomes('$DATA_PATH'); catch, exit(1); end; exit;"

echo "DONE"
