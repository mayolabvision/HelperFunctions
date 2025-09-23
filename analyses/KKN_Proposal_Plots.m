%% Proposal plot script

data_path = '/Users/kendranoneman/Data/sapu_dualhemi';

% get the filenames for all of the *.mat files in data_path
files = dir(fullfile(data_path, '*.mat'));

% initialize an empty cell array to store the contents of each .mat file
S_all = cell(1, numel(files));

% loop through each file
for i = 1:numel(files)
    fprintf('Loading file %d of %d: %s\n', i, numel(files), files(i).name);
    
    % load in each file, which contains a struct S
    tmp = load(fullfile(data_path, files(i).name), 'S');
    
    % add struct S to a cell array
    S_all{i} = tmp.S;
end

%%
