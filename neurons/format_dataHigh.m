function D = format_dataHigh(PATH, SESSION, varargin)
    % Converts Emily's data format to DataHigh-compatible format.
    %
    %%% Required inputs: %%%
    %   PATH    - Path to the data directory
    %             (e.g. PATH = '/Users/kendranoneman/Data/emily_data')
    %   SESSION - Session/file identifier w/out .mat file extension
    %             (e.g. SESSION = 's236')
    %
    %%% Optional parameters: %%%
    %   BLOCKS            - Array of blocks to process/concatenate into one struct, default is false (combines all blocks)
    %                       (e.g., BLOCKS = [1,2,3])
    %   CORRECT_ONLY      - Boolean, true to include only correct trials, false otherwise
    %                       (e.g., CORRECT_ONLY = true)
    %   CONDITION_BY      - String, 'block' or 'dir', default is 'block'
    %                       (e.g., CONDITION_BY = 'dir') 
    %   CONDITION_COLORS  - Cell array of colors, default is 12 different colors
    %                       (e.g., CONDITION_COLORS = {[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]})
    %   MIN_FR_HZ         - Minimum firing rate in Hz for neurons to include, default is 0
    %                       (e.g., MIN_FR_HZ = 1)
    %
    %%% Outputs: %%%
    %   D - Struct array compatible with DataHigh
    %
    %%% Example usage: %%%
    % D = format_dataHigh('/Users/kendranoneman/Data/emily_data', 's236', 'BLOCKS', [1,2,3], 'CORRECT_ONLY', true, ...
    %                     'CONDITION_BY', 'block', 'CONDITION_COLORS', {[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]}, ...
    %                     'MIN_FR_HZ', 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TO-DOS:
    % 1. Add 'alignCode' as input, for aligning spikes to either stimOnset,
    % saccadeOnset, etc...

    % Default values for optional parameters
    defaultBlocks = false; % Default to combine all blocks
    defaultCorrectOnly = false; % Default to include all trials
    defaultConditionBy = 'block'; % Default to condition by block
    defaultConditionColors = {[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], ...
                              [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75], ...
                              [1, 0.5, 0], [0.5, 1, 0], [0, 0.5, 1]}; % Default 12 colors
    defaultMinFrHz = 0; % Default minimum firing rate

    % Create an input parser
    p = inputParser;
    addRequired(p, 'PATH', @ischar); % PATH must be a string
    addRequired(p, 'SESSION', @ischar); % SESSION must be a string
    addParameter(p, 'BLOCKS', defaultBlocks, @(x) isnumeric(x) || (islogical(x) && ~x)); % BLOCKS must be numeric or false
    addParameter(p, 'CORRECT_ONLY', defaultCorrectOnly, @islogical); % CORRECT_ONLY must be logical
    addParameter(p, 'CONDITION_BY', defaultConditionBy, @ischar); % CONDITION_BY must be a string
    addParameter(p, 'CONDITION_COLORS', defaultConditionColors, @iscell); % CONDITION_COLORS must be a cell array
    addParameter(p, 'MIN_FR_HZ', defaultMinFrHz, @isnumeric); % MIN_FR_HZ must be numeric

    % Parse the inputs
    parse(p, PATH, SESSION, varargin{:});

    % Assign parsed values to variables
    PATH = p.Results.PATH;
    SESSION = p.Results.SESSION;
    BLOCKS = p.Results.BLOCKS;
    CORRECT_ONLY = p.Results.CORRECT_ONLY;
    CONDITION_BY = p.Results.CONDITION_BY;
    CONDITION_COLORS = p.Results.CONDITION_COLORS;
    MIN_FR_HZ = p.Results.MIN_FR_HZ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 1: Load data
    data = load('-mat', sprintf('%s/%s.mat', PATH, SESSION)); % Load the session data

    %% STEP 2: Find longest possible trial duration
    % DataHigh requires that all spike trains have the same length across all trials.
    % Therefore, we need to determine the longest possible trial first and use that length to bin all spike trains to ensure consistency.
    
    allTcodes  =  cellfun(@(x) vertcat(x.tcodes{:}), data.trialinfo, 'UniformOutput', false); % Concatenate tcodes from all trials
    allTcodes  =  vertcat(allTcodes{:}); % Concatenate all blocks
    edges      =  0:(1/1000):max(allTcodes(:,3)); % Define edges for histograms in 1ms bins

    %% STEP 3: Pull out fieldnames of "spikes_blocks"
    % In case datafiles in the future don't always have 6 spikes_blocks
    allFields  =  fieldnames(data.spikesbyclust(1)); % Get all field names from spikesbyclust
    spkBlocks  =  find(cellfun(@iscell, struct2cell(data.spikesbyclust(1))) == 1); % Find fields that are cell arrays
    blockNames =  allFields(spkBlocks); % Get names of those fields

    % If you don't define BLOCKS, this code will bin the spikes for each spikes_block separately and then concatenate them into one struct
    % But if you want to keep the spikes_blocks separate, then you need to specify that in the function input
    if ~isnumeric(BLOCKS) % If BLOCKS is not specified as numeric
        BLOCKS = 1:length(spkBlocks); % Use all blocks
    end
    blockNames = blockNames(BLOCKS); % Select block names based on BLOCKS

    %% STEP 4: Bin spikes for each spikes_block, cluster, and trial
    % You can specify if you only want to include correct trials or not 
    D_all = cell(1, length(BLOCKS));
    for block = 1:length(BLOCKS)
        D_eachBlock = struct([]); % Initialize empty struct for each block
        if CORRECT_ONLY % If only correct trials are included
            these_trials = find(logical(cellfun(@(q) sum(q(:,2)==150), data.trialinfo{block}.tcodes)) == 1); % Find correct trials
        else
            these_trials = 1:length(data.trialinfo{block}.tcodes); % Include all trials
        end
        for trial = 1:length(these_trials)
            spks_perTrial = zeros(length(data.spikesbyclust), length(edges) - 1); % Initialize spike train matrix for each trial
            for cluster = 1:length(data.spikesbyclust)
                spktimes = data.spikesbyclust(cluster).(blockNames{block}){these_trials(trial)}; % Get spike times for each cluster
                [spks_binned, ~] = histcounts(spktimes, edges); % Bin spike times into histogram
                spks_perTrial(cluster, :) = spks_binned; % Assign histogram to spike train matrix
            end

            D_eachBlock(trial).data = spks_perTrial; % Assign spike train matrix to data field
            if isequal(CONDITION_BY, 'block') % Condition by block
                conditionNum = block; % Condition number is block number
                D_eachBlock(trial).condition = blockNames{block}; % Condition name is block name
            elseif isequal(CONDITION_BY, 'dir') % Condition by direction
                thisTrial_dir = data.trialinfo{block}.dirs(these_trials(trial)); % Get direction for this trial
                conditionNum = find((sort(unique(data.trialinfo{block}.dirs)) == thisTrial_dir) == 1); % Find condition number for direction
                D_eachBlock(trial).condition = char(string(thisTrial_dir)); % Assign condition name as direction
            end
            D_eachBlock(trial).epochStarts = 1; % Start of epoch
            D_eachBlock(trial).epochColors = CONDITION_COLORS{conditionNum}; % Assign color for condition
        end
        D_all{block} = D_eachBlock; % Add block data to D_all
    end

    D = horzcat(D_all{:}); % Concatenate all blocks into one struct array

    % If you want to remove neurons that didn't fire enough, before saving the file
    % Trim the data field in each element of D using logical indexing
    numNeurons = size(D(1).data, 1);
    numTrials = length(D);
    meanFR_hz = zeros(numNeurons, 1);
    for trial = 1:numTrials
        meanFR_hz = meanFR_hz + sum(D(trial).data, 2) / (size(D(trial).data, 2) / 1000); % Sum spike counts for each neuron
    end
    meanFR_hz = meanFR_hz / numTrials; % Calculate mean firing rate
    
    neuronsToKeep = meanFR_hz > MIN_FR_HZ; % Use MIN_FR_HZ for threshold
    for i = 1:length(D)
        D(i).data = D(i).data(neuronsToKeep, :);
    end

    % Construct filename with details of included blocks, correct trials, condition type, and min firing rate
    blockStr = sprintf('blocks_%s', strjoin(string(BLOCKS), '_')); % Convert blocks to a string with underscores
    correctStr = 'allTrials'; % Default string for including all trials
    if CORRECT_ONLY
        correctStr = 'correctOnly'; % Change string if only correct trials are included
    end
    filename = sprintf('%s/dh_%s_%s_%s_%s_minFR%dHz.mat', PATH, SESSION, blockStr, correctStr, CONDITION_BY, MIN_FR_HZ); % Construct the full filename
    
    save(filename, 'D', '-v7.3'); % Save the struct array to a .mat file with the constructed filename
end