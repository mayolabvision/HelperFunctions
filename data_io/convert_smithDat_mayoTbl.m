function tbl = convert_smithDat_mayoTbl(dat,varargin)
    % format_dataTable - Processes neural and behavioral data from nev and ns5 files into a structured table.
    %
    % This function processes and formats data from neural and behavioral recording files. It takes in 
    % neural event (nev) and neural data (ns5) files, extracts relevant data based on task parameters, 
    % and returns a structured table with trial information, conditions, and eye movement data.
    %
    % The function supports optional parameters to customize the extraction of neural data, eye movement 
    % signals, and task-specific data. It can also format the output based on specific task types (e.g., 'rfmp', 
    % 'fstm', 'mdir', 'purs') and includes the option to organize data in a table format (because structs suck). 
    %
    %%%% Required inputs: %%%
    %   dat      -    Array of structs with the following fields:
    %       - block: The block number of the trial
    %       - time: Start and end times of the trial in seconds
    %       - text: Trial-specific text or annotations
    %       - trialcodes: Trial-specific codes for event markers
    %       - result: The result of the trial, such as "correct" or "incorrect"
    %       - params: The parameters associated with the trial (e.g., reaction times, crossing times)
    %       - eyes: Processed eye position data (in degrees)
    %       - pupil: Pupil data
    %       - diode: Diode signal data
    %       - spiketimes: Spike times for the specified neural channels, if applicable
    %       - net_labels: Spike sorting labels, if spike sorting is enabled
    %
    %%%% Optional parameters: %%%
    %   TASK_NAME   -  Name of the task, if you want to label it something other than what you named is in Ex ('sessionNumber').
    %
    %   LFP         -  Default = empty, if contains something should be cell array with same height as number of trials 
    %
    %%%% Outputs: %%%
    %   tbl  -  A table with many columns, including but not limited to:
    %           - trialName: Name or identifier of the trial
    %           - conditions: Conditions associated with the trial (e.g., jump, pursuitSpeed), with NaN if the condition is missing
    %           - result: Outcome of the trial, in categorical format (e.g., 'correct', 'incorrect')
    %           - params: Parameters of the trial and block (e.g., screen distance, target duration)
    %           - eyePos, eyeVel, eyeAcc: Filtered and smoothed eye movement data (position, velocity, and acceleration)
    %           - eyePos_raw: Raw, unfiltered eye data
    %           - spiketimes: Spike times for the specified neural channels, if applicable
    %           - net_labels: Spike sorting labels, if spike sorting is enabled
    %
    %%%% Example usage: %%%
    %   tbl = format_dataTable(nev, out_ns5, 'CONVERT_TO_TABLE', true)
    %
    % This will process the nev and ns5 files, with data organized in a table and task-specific data converted accordingly.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an input parser
    p = inputParser;
    addRequired(p, 'dat', @(x) (isnumeric(x)) || isstruct(x));
    addParameter(p, 'TASK_NAME', [], @ischar)
    addParameter(p, 'LFP', []);

    % Parse the inputs
    parse(p, dat, varargin{:});

    % Assign parsed values to variables
    dat = p.Results.dat;
    TASK_NAME = p.Results.TASK_NAME;
    LFP = p.Results.LFP;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(TASK_NAME)
        if isnumeric(dat(1).params.block.sessionNumber)
            TASK_NAME = sprintf('s%d',dat(1).params.block.sessionNumber);
        else
            TASK_NAME = dat(1).params.block.sessionNumber;
        end
    end

    tbl1 = struct2table(dat);
    tbl = table();

    % re-arranging table to be easier to access data 
    tbl.trialName = cellfun(@(q) [TASK_NAME, '.', sprintf('%04d', q)], num2cell(1:height(tbl1))', 'uni', 0);
    tbl.trialName = categorical(string(tbl.trialName));
    % tbl.block = tbl1.block;
    
    % Pull out conditions from trial names, separated by ';' delimeter
    if any(contains(TASK_NAME, {'rfmp', 'rfMapping'}))
        % Define the function to process each string
        process_string = @(input_str) [...
            str2double(cellfun(@(x) x{1}, regexp(input_str, 'xpos=([-0-9]+)', 'tokens'), 'UniformOutput', false))', ...
            str2double(cellfun(@(x) x{1}, regexp(input_str, 'ypos=([-0-9]+)', 'tokens'), 'UniformOutput', false))' ...
        ];

        % Apply the function to each cell in conditions
        conditions = cellfun(@(q) round(pix2deg(q, tbl1(1,:).params.block.screenDistance, tbl1(1,:).params.block.pixPerCM)), cellfun(process_string, tbl1.text, 'uni', 0), 'uni', 0);
        tbl.conditions = cellfun(@(q) num2cell(q,2), conditions, 'uni', 0);

    else
        pattern = '([^0-9;]+)(?==)';
        matches = cellfun(@(q) regexp(q, pattern, 'match'), tbl1.text, 'uni', 0);

        % Get all unique condition names from all the matches
        cols = sort(unique(horzcat(matches{:})));
        conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), tbl1.text, 'uni', 0);
        
        ordered_conditions = cell(size(conditions));
        for i = 1:length(conditions)
            temp = repmat({NaN}, 1, length(cols));
            
            current_row_cols = matches{i};
            for j = 1:length(current_row_cols)
                col_idx = find(strcmp(cols, current_row_cols{j}));
                
                temp{col_idx} = conditions{i}{j};
            end
            
            % Store the reordered conditions
            ordered_conditions{i} = temp;
        end
        
        % Convert the ordered conditions to a matrix
        conditions_matrix = cell2mat(vertcat(ordered_conditions{:}));
        
        % Create the table with the final ordered conditions
        for c = 1:length(cols)
            tbl.(cols{c}) = conditions_matrix(:, c);
        end
        % conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), tbl1.text, 'uni', 0);
        % conditions = cell2mat(vertcat(conditions{:}));
        % for c=1:length(cols)
        %     tbl.(cols{c}) = conditions(:,c);
        % end
    end

    if any(contains(TASK_NAME, {'purs','pursuit'}))
        tbl1.result(tbl1.result==0 | tbl1.result==154) = 167;
    end

    if iscell(tbl1.result)
        tbl.result = convertBetween_eventCodes_eventNames(tbl1.result);
    else
        tbl.result = convertBetween_eventCodes_eventNames(num2cell(tbl1.result));
    end
    tbl.result = categorical(string(tbl.result));

    % Make array of times of start/end time per trial, for aligning with trial codes and indexing eye data
    times_ms = cellfun(@(q) round(q(1)*1000:(q(2)+1)*1000), num2cell(tbl1.time,2), 'uni', 0);
    trialStarts = cellfun(@(q,r) find(q == round(r(r(:,2)==1,3)*1000)), times_ms, tbl1.trialcodes, 'uni', 0);
    eventCodes = tbl1.trialcodes; eventCodes = vertcat(eventCodes{:});
    eventCodes = num2cell(sort(unique(eventCodes(:,2))))';

    eventNames = convertBetween_eventCodes_eventNames(eventCodes);
    trialMarkers = cellfun(@(t) cellfun(@(q,r) find(ismember(q,round(r(r(:,2)==t,3)*1000))), times_ms, tbl1.trialcodes, 'uni', 0), eventCodes, 'uni', 0)';
    trialMarkers = horzcat(trialMarkers{:});
    trialMarkers(cellfun('isempty',trialMarkers)) = {NaN};

    for m=1:length(eventNames)
        if sum(cellfun(@(q) size(q,1)>1, trialMarkers(:,m), 'uni', 1))>0 || sum(cellfun(@(q) size(q,2)>1, trialMarkers(:,m), 'uni', 1))>0
            if numel(unique(cellfun(@numel ,trialMarkers(:,m),'uni',1))) == 1
                codeSplt = cellfun(@(q) num2cell(q), trialMarkers(:,m), 'uni', 0);
                tbl.(eventNames{m}) = vertcat(codeSplt{:});
            else
                tbl.(eventNames{m}) = trialMarkers(:,m);
            end
        else
            tbl.(eventNames{m}) = cell2mat(trialMarkers(:,m));
        end
    end

    tbl.params = tbl1.params;

    eyePos = cellfun(@(x,y) filterEyeTraces_EyeLink(x(:,y:end),'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',84,'PLOT_TRIAL',false), tbl1.eyedata, trialStarts, 'uni', 0);
    eyeVel = cellfun(@(q) calcDerivative_eyeTraces(q), cellfun(@(x,y) filterEyeTraces_EyeLink(x(:,y:end),'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',40,'PLOT_TRIAL',false), tbl1.eyedata, trialStarts, 'uni', 0), 'uni', 0);
    eyeAcc = cellfun(@(q) calcDerivative_eyeTraces(q), eyeVel, 'uni', 0);

    tbl.eyePos = eyePos; tbl.eyeVel = eyeVel; tbl.eyeAcc = eyeAcc;
    tbl.pupil = tbl1.pupil; tbl.diode = tbl1.diode; tbl.eyePos_raw = eyePos;

    if ismember('spiketimes', tbl1.Properties.VariableNames)
        tbl.spiketimes = tbl1.spiketimes; 
    end
    
    if ismember('net_labels', tbl1.Properties.VariableNames)
        tbl.net_labels = tbl1.net_labels;
    end

    if ~isempty(LFP)
        tbl.lfp = LFP;
    end

    %tbl.ns5_samps = tbl1.ns5_samps;


    %%%%%%%%%%%%%%%% TEMPORARY CODE FOR ANI %%%%%%%%%%%%%%%%%%%%%
    if ismember('SACCADE_BLOCK', tbl.Properties.VariableNames) && ismember('PURSUIT_BLOCK', tbl.Properties.VariableNames)
        tbl.blockType = cell(size(tbl, 1), 1);
        tbl.blockType(~isnan(tbl.SACCADE_BLOCK)) = {'mdir'};
        tbl.blockType(~isnan(tbl.PURSUIT_BLOCK)) = {'purs'};

        tbl = movevars(tbl,{'blockType'},'After','trialName');

        tbl.SACCADE_BLOCK = [];
        tbl.PURSUIT_BLOCK = [];

        tbl.blockType = categorical(string(tbl.blockType));
        tbl.distance(tbl.blockType=='pursuit') = NaN;
        tbl.jump(tbl.blockType=='saccade') = NaN;
        tbl.pursuitSpeed(tbl.blockType=='saccade') = NaN;
        tbl.pursuit_fixDuration(tbl.blockType=='saccade') = NaN;

        % changing distance to be in deg
        tbl.distance = cellfun(@(q) round(pix2deg(q, tbl(1,:).params.block.screenDistance, tbl(1,:).params.block.pixPerCM)), num2cell(tbl.distance), 'uni', 1);

        % Create empty columns for pursuit-related variables
        tbl.pursuitOnset = nan(height(tbl), 1);
        tbl.pursuitLatency = nan(height(tbl), 1);
        tbl.msOffset = nan(height(tbl), 1);
        tbl.CROSSING_TIME = nan(height(tbl), 1);
        tbl.csTimes = nan(height(tbl), 3);
        tbl.csVelocity = nan(height(tbl), 1);
        tbl.csAngle = nan(height(tbl), 1);
        tbl.pursType = repmat({'NaN'}, height(tbl), 1);

        rowsForPursuit = find(tbl.blockType=='purs');
        
        % Loop over only the relevant rows for pursuit
        for t = find(rowsForPursuit)'
            if isequal(tbl.result(t), "CORRECT") && length(tbl.TARG_ON{t})==1
                [pursuit_onset, rxnTime, msOffset, csOnset, csVelocity, csPeak, csOffset, csAngle, csType] = detect_pursuitOnset(tbl.eyePos{t}, tbl.eyeVel{t}, tbl.TARG_ON{t}, tbl(t,:).params.block.crossingTime, tbl.pursuitSpeed(t), tbl.newAngle(t), 'PLOT_TRACES', false);
                
                tbl.pursuitOnset(t) = pursuit_onset;
                tbl.pursuitLatency(t) = rxnTime;
                tbl.msOffset(t) = msOffset;
                tbl.csTimes(t, :) = [csOnset, csPeak, csOffset];
                tbl.csVelocity(t) = csVelocity;
                tbl.csAngle(t) = csAngle;
                tbl.pursType{t} = csType;
            else
                tbl.pursType{t} = 'NaN';
            end
            
            if isequal(tbl.result(t), 'CORRECT')
                tbl.CROSSING_TIME(t) = tbl(t,:).params.block.crossingTime;
            end
        end
        
        tbl.pursType = categorical(string(tbl.pursType));
        
        tbl = movevars(tbl,{'pursuitOnset','pursuitLatency','msOffset','pursType','csTimes','csVelocity','csAngle'},'Before','result');
        tbl = movevars(tbl,{'CROSSING_TIME'},'After','TARG_ON');
        
    end

    %%%%%%%%%%%%% task-specific re-arranging and calculations %%%%%%%%%%%%%
    if any(contains(TASK_NAME, {'mdir', 'dirmem'}))
        tbl(:, 'MEM_GUIDED_SACC') = [];

        if ~ismember('distance', tbl.Properties.VariableNames)
            tbl.distance = cellfun(@(q) q.distance, {tbl.params.block}.', 'uni', 1);
            tbl = movevars(tbl,{'distance'},'After','angle');
        end
        tbl.distance = cellfun(@(q) round(pix2deg(q,tbl(1,:).params.block.screenDistance,tbl(1,:).params.block.pixPerCM)), num2cell(tbl.distance), 'uni', 1);

        if all(ismember({'targetOnsetDelay', 'delay'}, tbl.Properties.VariableNames))
            tbl.fixDuration = tbl.targetOnsetDelay+tbl(1,:).params.block.targetDuration+tbl.delay;

            tbl = movevars(tbl,{'targetOnsetDelay','delay','fixDuration'},'Before','result');
        end

        if ~iscell(tbl.FIX_OFF) 
            tbl.FIX_OFF = num2cell(tbl.FIX_OFF);
        end

        tbl.saccLatency = tbl.SACCADE-cellfun(@(q) q(1), tbl.FIX_OFF);
             

    elseif any(contains(TASK_NAME, {'purs','pursuit'}))
        % Define the columns to replace and their new names
        cols_to_replace = {'TARG_ON', 'PURSUIT_TARG', 'TARG_OFF'};
        new_names = {'PURSUIT_TARG_ON', 'PURSUIT_TARG_ON', 'PURSUIT_TARG_OFF'};
        
        % Loop through each column to check and replace
        for i = 1:numel(cols_to_replace)
            if ismember(cols_to_replace{i}, tbl.Properties.VariableNames)
                tbl.Properties.VariableNames{ismember(tbl.Properties.VariableNames, cols_to_replace{i})} = new_names{i};
            end
        end

        if ~ismember('pursuitSpeed',tbl.Properties.VariableNames)
            tbl.pursuitSpeed = repmat(tbl(1,:).params.block.pursuitSpeed,height(tbl),1);
            tbl = movevars(tbl,{'pursuitSpeed'},'Before','fixDuration');
        end
        if ~ismember('jump',tbl.Properties.VariableNames)
            if isfield(tbl(1,:).params.block,'jump')
                tbl.jump = repmat(tbl(1,:).params.block.jump,height(tbl),1);
            else
                tbl.jump = zeros(height(tbl),1);
            end
            tbl = movevars(tbl,{'jump'},'Before','fixDuration');
        end

        [pursuitOnsets,rxnTimes,msOffsets,csOnsets,csVelocities,csPeaks,csOffsets,csAngles,crossingTimes] = deal(nan(height(tbl), 1));
        csTypes = cell(height(tbl),1);
        for t = 1:height(tbl)
            if isequal(tbl.result(t),"CORRECT")
                
                if isfield(tbl(1,:).params.block,'crossingTime')
                    [pursuit_onset,rxnTime,msOffset,csOnset,csVelocity,csPeak,csOffset,csAngle,csType] = detect_pursuitOnset(tbl.eyePos{t},tbl.eyeVel{t},tbl.PURSUIT_TARG_ON(t),tbl(t,:).params.block.crossingTime,tbl.pursuitSpeed(t),tbl.angle(t),'PLOT_TRACES',false);
                else
                    [pursuit_onset,rxnTime,msOffset,csOnset,csVelocity,csPeak,csOffset,csAngle,csType] = detect_pursuitOnset(tbl.eyePos{t},tbl.eyeVel{t},tbl.PURSUIT_TARG_ON(t),110,tbl.pursuitSpeed(t),tbl.angle(t),'PLOT_TRACES',false);
                end
                pursuitOnsets(t) = pursuit_onset; rxnTimes(t) = rxnTime; msOffsets(t) = msOffset; csOnsets(t) = csOnset; csVelocities(t) = csVelocity; csPeaks(t) = csPeak; csOffsets(t) = csOffset; csAngles(t) = csAngle; csTypes{t} = csType;
            else
                csTypes{t} = 'NaN';
            end

        end

        tbl.pursuitOnset = pursuitOnsets; tbl.pursuitLatency = rxnTimes;
        tbl.msOffset = msOffsets; tbl.CROSSING_TIME = crossingTimes;
        tbl.csTimes = [csOnsets, csPeaks, csOffsets]; tbl.csVelocity = csVelocities; tbl.csAngle = csAngles;
        tbl.pursType = csTypes; tbl.pursType = categorical(string(tbl.pursType));
        
        tbl = movevars(tbl,{'pursuitOnset','pursuitLatency','msOffset','pursType','csTimes','csVelocity','csAngle'},'Before','result');
        tbl = movevars(tbl,{'CROSSING_TIME'},'After','PURSUIT_TARG_ON');

    % elseif any(contains(TASK_NAME, {'rfmp','rfMapping'})) && ismember('STIM_ON', tbl.Properties.VariableNames)
    %     tbl.STIM_ON(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_ON(tbl.result~='CORRECT'), 'uni', 0);
    %     tbl.conditions(tbl.result~="CORRECT") = cellfun(@(q) q(1:end-1), tbl.conditions(tbl.result~='CORRECT'), 'uni', 0);
    % end

    if ismember('emptyCnd', tbl.Properties.VariableNames)
        tbl.emptyCnd = [];
    end

    if ismember('STIM8_ON', tbl.Properties.VariableNames)
        tbl.STIM8_ON = [];
    end

    if ismember('ALIGN_PULSE', tbl.Properties.VariableNames)
        tbl = movevars(tbl,{'ALIGN_PULSE'},'After','START_TRIAL');
    end

    if ismember('IGNORED', tbl.Properties.VariableNames) && ismember('CORRECT', tbl.Properties.VariableNames)
        tbl = movevars(tbl,{'IGNORED'},'After','CORRECT');
    end

    % if any(contains(TASK_NAME, {'rfmp', 'rfMapping'}))
    %     tbl = tbl(~cellfun(@(q) any(isnan(q)), tbl.STIM_OFF, 'uni', 1),:);
    % end

end
