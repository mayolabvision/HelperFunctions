#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=mtlab
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --error=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-02:00:00

echo "Job started at $(date)"

module purge
module load matlab/R2023a
conda activate /ihome/pmayo/knoneman/.conda/envs/npy2mat

# Assign variables from positional parameters
SESSION="$1"
RUN_TYPE="${2:-unleashed}"
SWEEP_NAME="${3:-none}"

echo "SESSION: $SESSION"
echo "RUN_TYPE: $RUN_TYPE"

KILO4_PATH0="/ix1/pmayo/lab_NHPdata/${SESSION}/${SESSION}_imec0/kilosort4_${RUN_TYPE}"
KILO4_PATH1="/ix1/pmayo/lab_NHPdata/${SESSION}/${SESSION}_imec1/kilosort4_${RUN_TYPE}"

if [[ "$RUN_TYPE" == "sweep" ]]; then
    KILO4_PATH0="${KILO4_PATH0}/${SWEEP_NAME}"
    KILO4_PATH1="${KILO4_PATH1}/${SWEEP_NAME}"
fi
echo "KILO4_PATH0: $KILO4_PATH0"


# Define the varargin parameters
RAW_PATH='/ix1/pmayo/lab_NHPdata'
OUT_PATH='/ix1/pmayo/lab_NHPdata'
NEV_PATH='/ihome/pmayo/knoneman/Packages/nevutils'
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'
PROBE_TYPE='np'
PARSE_KS="true"

echo "PROBE_TYPE: $PROBE_TYPE"
echo "PARSE_KS: $PARSE_KS"
echo "RUN_TYPE: $RUN_TYPE"

##########################################################################3

# Convert .npy files in kilosort4 directory to .mat
python -c "
import sys
sys.path.append('/ihome/pmayo/knoneman/Packages/HelperFunctions/utils')
from convert_npy_to_mat import convert_npy_to_mat
print('Running convert_npy_to_mat on: $KILO4_PATH0')
convert_npy_to_mat('$KILO4_PATH0')
"

python -c "
import sys
import os
sys.path.append('/ihome/pmayo/knoneman/Packages/HelperFunctions/utils')
from convert_npy_to_mat import convert_npy_to_mat

kilo_path = '$KILO4_PATH1'
if os.path.exists(kilo_path):
    print(f'Running convert_npy_to_mat on: {kilo_path}')
    convert_npy_to_mat(kilo_path)
else:
    print(f'Skipping convert_npy_to_mat: {kilo_path} does not exist.')
"

# First MATLAB call: run process_fullRecording
# First MATLAB call: run process_fullRecording
matlab -nodisplay <<EOF
try
    addpath(genpath('$HELPERS_PATH'));
    fprintf('Running process_fullRecording for $1\n');
    process_fullRecording('${SESSION}', ...
        'RAW_DATA_PATH', '$RAW_PATH', ...
        'OUT_DATA_PATH', '$OUT_PATH', ...
        'NEVUTIL_PATH', '$NEV_PATH', ...
        'PROBE_TYPE', '$PROBE_TYPE', ...
        'PARSE_KILOSORT', evalin('base', '$PARSE_KS'), ...
        'RUN_TYPE', '$RUN_TYPE', ...
        'SWEEP_NAME', '$SWEEP_NAME');
catch err
    disp('ERROR in process_fullRecording:');
    disp(getReport(err));
    exit(1);
end
exit
EOF

#matlab -nodisplay -r "try, addpath(genpath('$HELPERS_PATH')); ia_trialOutcomes('$DATA_PATH'); catch, exit(1); end; exit;"

echo "DONE"
