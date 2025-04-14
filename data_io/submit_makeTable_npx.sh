#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --job-name=mtlab
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --error=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-00:09:59

echo "Job started at $(date)"

module purge
module load matlab/R2023a
conda activate /ihome/pmayo/knoneman/.conda/envs/npy2mat

# Define the varargin parameters
RAW_PATH='/ix1/pmayo/lab_NHPdata'
OUT_PATH='/ix1/pmayo/lab_NHPdata'
CSV_PATH='/ix1/pmayo/lab_NHPdata/RECORDING_INFO.csv'
NEV_PATH='/ihome/pmayo/knoneman/Packages/nevutils'
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'
KILO4_PATH0="/ix1/pmayo/lab_NHPdata/${1}/${1}_imec0/kilosort4"
KILO4_PATH1="/ix1/pmayo/lab_NHPdata/${1}/${1}_imec1/kilosort4"
DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}.mat"

echo "Using DATA_PATH: $DATA_PATH"

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
sys.path.append('/ihome/pmayo/knoneman/Packages/HelperFunctions/utils')
from convert_npy_to_mat import convert_npy_to_mat
print('Running convert_npy_to_mat on: $KILO4_PATH1')
convert_npy_to_mat('$KILO4_PATH1')
"

# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
try
    addpath(genpath('$HELPERS_PATH'));
    fprintf('Running process_NeuropixRecording_KKN for $1\n');
    process_NeuropixRecording_KKN('$DATA_PATH', ...
        'RAW_DATA_PATH', '$RAW_PATH', ...
        'OUT_DATA_PATH', '$OUT_PATH', ...
        'RECD_CSV_PATH', '$CSV_PATH', ...
        'NEVUTIL_PATH', '$NEV_PATH');
catch err
    disp('ERROR in process_NeuropixRecording_KKN:');
    disp(getReport(err));
    exit(1);
end
exit
EOF

#matlab -nodisplay -r "try, process_NeuropixRecording_KKN('$1', '$2', '$3', 'RAW_DATA_PATH', '$RAW_PATH', 'OUT_DATA_PATH', '$OUT_PATH', 'RECD_CSV_PATH', '$CSV_PATH', 'NASNET_PATH', '$NET_PATH', 'NEVUTIL_PATH', '$NEV_PATH'); catch, exit(1); end; exit;"

#matlab -nodisplay -r "try, ia_trialOutcomes('$DATA_PATH'); catch, exit(1); end; exit;"

echo "DONE"
