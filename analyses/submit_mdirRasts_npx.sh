#!/bin/bash -l
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --job-name=mdir
#SBATCH --error=/ix1/pmayo/matlab/outfiles/error_%A.out
#SBATCH --output=/ix1/pmayo/matlab/outfiles/out_%A.out
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

# Define the varargin parameters
HELPERS_PATH='/ihome/pmayo/knoneman/Packages/HelperFunctions'
DATA_PATH="/ix1/pmayo/lab_NHPdata/${1}/${1}.mat"

# First MATLAB call: run process_NeuropixRecording_KKN
matlab -nodisplay <<EOF
try
    addpath(genpath('$HELPERS_PATH'));
    fprintf('Running ia_mdirRasters for $1\n');
    ia_mdirRasters('$DATA_PATH', ...
        'IMEC', ${2}, ...
        'ALIGN', '${3}');
catch err
    disp('ERROR in ia_mdirRasters:');
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
