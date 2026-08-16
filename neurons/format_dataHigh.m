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
    %   ALIGN_BY          - Trial code (defined in exGlobals) to align the spike trains to
    %                       (e.g., ALIGN_BY = 'SACCADE' or ALIGN_BY = 255)
    %   SPK_INTERVAL      - Range around the ALIGN_BY code to pull out spikes, default is [-500, 500]
    %                       (e.g., SPK_INTERVAL = [-500, 500])
    %   BLOCKS            - Array of blocks to process/concatenate into one struct, default is false (combines all blocks)
    %                       (e.g., BLOCKS = [1,2,3])
    %   CORRECT_ONLY      - Boolean, true to include only correct trials, false otherwise
    %                       (e.g., CORRECT_ONLY = true)
    %   CONDITION_BY      - String, 'block' or 'dir', default is 'block'
    %                       (e.g., CONDITION_BY = 'dir') 
    %   CONDITION_GROUPS  - Cell array of block numbers or directions you wish to group together, default is false (none are grouped)
    %                       (e.g., CONDITION_GROUPS = {[1,3,5],[2,4,6]})
    %   CONDITION_NAMES   - Cell array of condition names
    %                       (e.g., CONDITION_NAMES = {'c1','c2'})
    %   CONDITION_COLORS  - Cell array of colors, default is 2 colors per condition 
    %                        (e.g., CONDITION_COLORS = {[[102,178,255]./255;[0,102,204]./255], [[255,102,102]./255;[204,0,0]./255]})
    %   EPOCH_STARTS      - Array of points to put markers 
    %                        (e.g., [-200,-100,0,100,200])
    %
    %%% Outputs: %%%
    %   D - Struct array compatible with DataHigh
    %
    %%% Example usage: %%%
    %   D = format_dataHigh('/Users/kendranoneman/Data/emily_data', 's236', 'ALIGN_BY', 'SACCADE', 'SPK_INTERVAL', [-500, 500], ...
    %                     'BLOCKS', [1,2,3], 'CORRECT_ONLY', true, 'CONDITION_BY', 'block', 'CONDITION_GROUPS', {[1,3,5],[2,4,6]}, 'CONDITION_NAMES', {'c1','c2'}...
    %                     'CONDITION_COLORS', {[[102,178,255]./255; [0,102,204]./255], [[255,102,102]./255; [204,0,0]./255]]}, ...
    %                     'EPOCH_STARTS',[-200,-100,0,100,200]);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addpath(genpath(fullfile(pwd, '..'))); % so you have access to the convertBetween_eventCodes_eventNames function

    % Default values for optional parameters
    defaultBlocks = false; % Default to combine all blocks
    defaultCorrectOnly = true; % Default to include all trials
    defaultAlignBy = 'SACCADE'; % Trial code or code name to align to
    defaultSpkInterval = [-300, 200]; % Default spike interval in ms
    defaultConditionBy = 'block'; % Default to condition by block
    defaultConditionGroups = false; % Default is to not group any of the conditions together
    defaultConditionNames = false; % Default is to name conditions based on existing values
    defaultConditionColors = repmat({[[102, 178, 255]./255;[0, 102, 204]./255], [[255, 102, 102]./255;[204, 0, 0]./255]}, 1, 6);
    defaultEpochStarts = false; % Default is to not specify any epochStarts
    

    % Create an input parser
    p = inputParser;
    addRequired(p, 'PATH', @ischar); % PATH must be a string
    addRequired(p, 'SESSION', @ischar); % SESSION must be a string
    addParameter(p, 'ALIGN_BY', defaultAlignBy, @(x) ischar(x) || isstring(x) || isnumeric(x)); % ALIGN_BY must be a string/char or double
    addParameter(p, 'SPK_INTERVAL', defaultSpkInterval, @(x) isnumeric(x) && numel(x) == 2); % SPK_INTERVAL must be a numeric array of size 2
    addParameter(p, 'BLOCKS', defaultBlocks, @(x) isnumeric(x) || (islogical(x) && ~x)); % BLOCKS must be numeric or false
    addParameter(p, 'CORRECT_ONLY', defaultCorrectOnly, @islogical); % CORRECT_ONLY must be logical
    addParameter(p, 'CONDITION_BY', defaultConditionBy, @ischar); % CONDITION_BY must be a string
    addParameter(p, 'CONDITION_GROUPS', defaultConditionGroups, @iscell); % CONDITION_GROUPS must be a cell array
    addParameter(p, 'CONDITION_NAMES', defaultConditionNames, @iscell); % CONDITION_GROUPS must be a cell array
    addParameter(p, 'CONDITION_COLORS', defaultConditionColors, @iscell); % CONDITION_COLORS must be a cell array
    addParameter(p, 'EPOCH_STARTS', defaultEpochStarts, @(x) isnumeric(x) || (islogical(x) && ~x)); % BLOCKS must be numeric or false

    % Parse the inputs
    parse(p, PATH, SESSION, varargin{:});

    % Assign parsed values to variables
    PATH = p.Results.PATH;
    SESSION = p.Results.SESSION;
    ALIGN_BY = p.Results.ALIGN_BY;
    SPK_INTERVAL = p.Results.SPK_INTERVAL;
    BLOCKS = p.Results.BLOCKS;
    CORRECT_ONLY = p.Results.CORRECT_ONLY;
    CONDITION_BY = p.Results.CONDITION_BY;
    CONDITION_GROUPS = p.Results.CONDITION_GROUPS;
    CONDITION_NAMES = p.Results.CONDITION_NAMES;
    CONDITION_COLORS = p.Results.CONDITION_COLORS;
    EPOCH_STARTS = p.Results.EPOCH_STARTS;

    % Check SPK_INTERVAL validity
    if SPK_INTERVAL(1) >= SPK_INTERVAL(2)
        error('SPK_INTERVAL:InvalidRange', 'The first number in SPK_INTERVAL must be less than the second.');
    end

    % Check that CONDITION_NAMES and CONDITION_GROUPS are either both false or have the same dimensions
    if ~(isequal(CONDITION_NAMES, false) && isequal(CONDITION_GROUPS, false)) && ...
       ~isequal(size(CONDITION_NAMES), size(CONDITION_GROUPS))
        error('CONDITION_NAMES and CONDITION_GROUPS must either both be false or have the same dimensions.');
    end

    % If you only have one color defined per condition, then D(i).epochStarts = 1
    if islogical(EPOCH_STARTS) && isequal(EPOCH_STARTS, false)
        if size(CONDITION_COLORS{1},1) == 1
            epoch_start = 1;
        else
            epoch_start = [1,abs(SPK_INTERVAL(1))];
        end

    elseif isnumeric(EPOCH_STARTS)
        EPOCH_STARTS = EPOCH_STARTS(EPOCH_STARTS<SPK_INTERVAL(2) & EPOCH_STARTS>=SPK_INTERVAL(1));
        % Sort, remove leading 0, and add 1 if necessary
        epoch_start = sort(EPOCH_STARTS + abs(SPK_INTERVAL(1)));
        epoch_start = epoch_start(epoch_start ~= 0);
        if epoch_start(1) ~= 1, epoch_start = [1, epoch_start]; end

        % Function to handle repeating colors
        repeatColors = @(colors, len) repmat(colors, ceil(len / size(colors, 1)), 1);
        
        % Loop through each cell in CONDITION_COLORS
        for i = 1:length(CONDITION_COLORS)
            colors = CONDITION_COLORS{i};
            CONDITION_COLORS{i} = repeatColors(colors, length(epoch_start)); % Repeat the colors
            CONDITION_COLORS{i} = CONDITION_COLORS{i}(1:length(epoch_start), :); % Trim to ensure exact length
        end
    else
        error('EPOCH_STARTS must be either false or an array of doubles.');
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 1: Load data
    data = load('-mat', sprintf('%s/%s.mat', PATH, SESSION)); % Load the session data

    %% STEP 2: Define duration of spike trains, based on SPK_INTERVAL
    % DataHigh requires that all spike trains have the same length across all trials.
    % Therefore, we need to determine the longest possible trial first and use that length to bin all spike trains to ensure consistency.
    
    allTcodes  =  cellfun(@(x) vertcat(x.tcodes{:}), data.trialinfo, 'UniformOutput', false); % Concatenate tcodes from all trials
    allTcodes  =  vertcat(allTcodes{:}); % Concatenate all blocks
    edges      =  (SPK_INTERVAL(1)/1000):(1/1000):(SPK_INTERVAL(2)/1000); % Define edges for histograms in 1ms bins

    %% STEP 3: Define which trial code to align spike trains to
    % The ALIGN_BY input can be either the word for the code (e.g., 'SACCADE') or the numeric code (e.g., 255)
    trialCodes = sort(unique(allTcodes(:,2)));
    possibleCodes = trialCodes(trialCodes>0 & trialCodes<1000);
    possibleCodeNames = convertBetween_eventCodes_eventNames(num2cell(possibleCodes));
    if isequal(class(ALIGN_BY), 'char') || isequal(class(ALIGN_BY), 'string')
        if ismember(ALIGN_BY, possibleCodeNames)
            alignCode = possibleCodes(ismember(possibleCodeNames, ALIGN_BY));
            alignName = ALIGN_BY;
        else
            fprintf('Error: ALIGN_BY "%s" does not match any possible code names.\n', ALIGN_BY);
            error('ALIGN_BY does not match any possible code names.');
        end
    elseif isequal(class(ALIGN_BY), 'double')
        if ismember(ALIGN_BY, possibleCodes)
            alignCode = ALIGN_BY;
            alignName = possibleCodeNames(possibleCodes==ALIGN_BY);
            alignName = alignName{1};
        else
            fprintf('Error: ALIGN_BY "%d" does not match any possible codes.\n', ALIGN_BY);
            error('ALIGN_BY does not match any possible codes.');
        end
    else
        fprintf('Error: ALIGN_BY is of unsupported type "%s".\n', class(ALIGN_BY));
        error('ALIGN_BY is of unsupported type.');
    end

    %% STEP 4: Pull out fieldnames of "spikes_blocks"
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

    %% STEP 5: Bin spikes for each spikes_block, cluster, and trial
    % You can specify if you only want to include correct trials or not in the function inputs
    % If the code you want to align spikes to doesn't exist for a given trial, it skips that trial
    D_all = cell(1, length(BLOCKS));
    for block = 1:length(BLOCKS)
        D_eachBlock = struct([]); % Initialize empty struct for each block
        if CORRECT_ONLY % If only correct trials are included
            these_trials = find(cellfun(@(q) sum(q(:,2)==150) + sum(q(:,2)==alignCode), data.trialinfo{block}.tcodes) == 2); % Find correct trials
        else
            these_trials = find(cellfun(@(q) sum(q(:,2)==alignCode), data.trialinfo{block}.tcodes) == 1);
        end
        %these_trials = 
        for trial = 1:length(these_trials)
            alignTime = data.trialinfo{block}.tcodes{these_trials(trial)}(data.trialinfo{block}.tcodes{these_trials(trial)}(:,2)==alignCode,3);
            spks_perTrial = zeros(length(data.spikesbyclust), length(edges) - 1); % Initialize spike train matrix for each trial
            for cluster = 1:length(data.spikesbyclust)
                spktimes = data.spikesbyclust(cluster).(blockNames{block}){these_trials(trial)}; % Get spike times for each cluster
                if ~isempty(spktimes)
                    [spks_binned, ~] = histcounts(spktimes-alignTime, edges); % Bin spike times into histogram
                else
                    spks_binned = zeros(1,size(edges,2)-1);
                end
                spks_perTrial(cluster, :) = spks_binned; % Assign histogram to spike train matrix
            end

            D_eachBlock(trial).data = spks_perTrial; % Assign spike train matrix to data field

            if isequal(CONDITION_BY, 'block') % Condition by block
                conditionNum = block; % Condition number is block number
                if iscell(CONDITION_NAMES)
                    D_eachBlock(trial).condition = CONDITION_NAMES{cellfun(@(q) ismember(conditionNum,q), CONDITION_GROUPS)};
                    D_eachBlock(trial).epochColors = CONDITION_COLORS{cellfun(@(q) ismember(conditionNum,q), CONDITION_GROUPS)};
                elseif isequal(CONDITION_NAMES, false)
                    D_eachBlock(trial).condition = blockNames{block}; % Condition name is block name
                    D_eachBlock(trial).epochColors = CONDITION_COLORS{conditionNum};
                else
                    error('CONDITION_NAMES must be either false or a cell array.');
                end
            elseif isequal(CONDITION_BY, 'dir') % Condition by direction
                thisTrial_dir = data.trialinfo{block}.dirs(these_trials(trial)); % Get direction for this trial
                conditionNum = sort(unique(data.trialinfo{block}.dirs)) == thisTrial_dir; % Find condition number for direction
                if iscell(CONDITION_NAMES)
                    D_eachBlock(trial).condition = CONDITION_NAMES{cellfun(@(q) ismember(thisTrial_dir,q), CONDITION_GROUPS)};
                    D_eachBlock(trial).epochColors = CONDITION_COLORS{cellfun(@(q) ismember(thisTrial_dir,q), CONDITION_GROUPS)};
                elseif isequal(CONDITION_NAMES, false)
                    D_eachBlock(trial).condition = char(string(thisTrial_dir)); % Assign condition name as direction
                    D_eachBlock(trial).epochColors = CONDITION_COLORS{conditionNum};
                else
                    error('CONDITION_NAMES must be either false or a cell array.');
                end
            end
            D_eachBlock(trial).type = 'traj';
            D_eachBlock(trial).epochStarts = epoch_start;
        end
        D_all{block} = D_eachBlock; % Add block data to D_all
    end

    D = horzcat(D_all{:}); % Concatenate all blocks into one struct array

    %% STEP 6: Save struct, D, to a .mat file
    % Construct filename with details of included blocks, correct trials, condition type, min firing rate, SPK_INTERVAL, and ALIGN_BY
    blockStr = sprintf('blocks_%s', strjoin(string(BLOCKS), '_')); % Convert blocks to a string with underscores
    correctStr = 'allTrials'; % Default string for including all trials
    if CORRECT_ONLY
        correctStr = 'correctOnly'; % Change string if only correct trials are included
    end
    spkIntervalStr = sprintf('spkInt_%d_%dms', SPK_INTERVAL(1), SPK_INTERVAL(2)); % Convert SPK_INTERVAL to a string
    alignByStr = sprintf('alignBy_%s', alignName); % Convert ALIGN_BY to a string
    filename = sprintf('%s/dh_%s_%s_%s_%s_%s_%s.mat', PATH, SESSION, blockStr, alignByStr, spkIntervalStr, correctStr, CONDITION_BY); % Construct the full filename
    
    save(filename, 'D', '-v7.3'); % Save the struct array to a .mat file with the constructed filename
end