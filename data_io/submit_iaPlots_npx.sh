#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --job-name=plts
#SBATCH --error=/ix1/pmayo/matlab/outfiles/error_%A.err
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=fail
#SBATCH --mail-user=knoneman@pitt.edu
#SBATCH --time=0-00:09:59

echo "Job started at $(date)"

module purge
module load matlab/R2023a

# Define the varargin parameters
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'
DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}.mat"

# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
try
    addpath(genpath('$HELPERS_PATH'));
    fprintf('Running ia_rfMaps for $1\n');
    ia_rfMaps('${1}', ...
        'IMEC', ${2});
catch err
    disp('ERROR in ia_rfMaps:');
    disp(getReport(err));
    exit(1);
end
exit
EOF

#if [ "$2" == "rfmp" ]; then
#    matlab -nodisplay -r "try, addpath('$HELPERS_PATH'); ia_rfMaps('$DATA_PATH', 'IMEC', '$3'); catch, exit(1); end; exit;"
#elif [ "$2" == "mdir" ]; then
#    matlab -nodisplay -r "try, addpath('$HELPERS_PATH'); ia_mdirRasters('$DATA_PATH', 'IMEC', '$3', 'ALIGN', '$4'); catch, exit(1); end; exit;"
#else
#    echo "Unknown function label: $3"
#    exit 1
#fi

echo "DONE"
