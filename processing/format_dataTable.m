function tbl = format_dataTable(nev,out,task_name,varargin)
    % blah
    %
    %%% Required inputs: %%%
    %   nev    -    nev file
    %   out    -    out file
    %
    %%% Optional parameters: %%%
    %   EYE_CHAN_LABELS  -  {Eye_HE, Eye_VE, DIODE, PUPIL}
    %                      (Default: {'10241', '10242', '10243', '10244'})
    
    %%% Outputs: %%%
    %   tbl  -  3-column array containing all spikes and trial codes
    %           (nev(:,1) = channel number, nev(:,2) = sort code, nev(:,3) = time w/ 30 kHz Fs)
    %            if nev(t,1) == 0 then nev(t,2) contains trial code 
    %
    %%% Example usage: %%%
    %   tbl = format_dataTable(nev,out,'EYE_CHAN_LABELS',{'10241', '10242', '10243', '10244'},'CONVERT_TO_TABLE',true)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Default values for optional parameters
    defaultEyeChanLabels = {'10241', '10242', '10243', '10244'};
    

    % Create an input parser
    p = inputParser;
    addRequired(p, 'nev', @(x) (isnumeric(x) && size(x, 2) == 3) || isstruct(x));
    addRequired(p, 'out', @isstruct);
    addParameter(p, 'EYE_CHAN_LABELS', defaultEyeChanLabels, (@(x) iscell(x) && length(x)==4)); % channel labels

    % Parse the inputs
    parse(p, nev, out, varargin{:});

    % Assign parsed values to variables
    nev = p.Results.nev;
    out = p.Results.out;
    EYE_CHAN_LABELS = p.Results.EYE_CHAN_LABELS;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dat = format_datTrials(nev,out,EYE_CHAN_LABELS);

    tbl1 = struct2table(dat);
    tbl = table();

    % re-arranging table to be easier to access data 
    tbl.trialName = cellfun(@(q) [task_name,char('.'),char(string(q))], num2cell(1:height(tbl1))','uni',0);
    tbl.trialName = categorical(string(tbl.trialName));
    tbl.block = tbl1.block;
    
    % Pull out conditions from trial names, separated by ';' delimeter
    if contains(task_name, 'rfmp')
        % Define the function to process each string
        process_string = @(input_str) [...
            str2double(cellfun(@(x) x{1}, regexp(input_str, 'xpos=([-0-9]+)', 'tokens'), 'UniformOutput', false))', ...
            str2double(cellfun(@(x) x{1}, regexp(input_str, 'ypos=([-0-9]+)', 'tokens'), 'UniformOutput', false))' ...
        ];

        % Apply the function to each cell in conditions
        conditions = cellfun(process_string, tbl1.text, 'uni', 0);
        tbl.conditions = cellfun(@(q) num2cell(q,2), conditions, 'uni', 0);
    else
        pattern = '([^0-9;]+)(?==)';
        matches = cellfun(@(q) regexp(q, pattern, 'match'), tbl1.text, 'uni', 0);
        cols = matches{1};
        conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), tbl1.text, 'uni', 0);
        conditions = cell2mat(vertcat(conditions{:}));
        for c=1:length(cols)
            tbl.(cols{c}) = conditions(:,c);
        end
    end

    if contains(task_name, 'purs')
        tbl1.result(tbl1.result==0 | tbl1.result==154) = 167;
    end

    tbl.result = convertBetween_eventCodes_eventNames(num2cell(tbl1.result));
    tbl.result = categorical(string(tbl.result));

    % Make array of times of start/end time per trial, for aligning with trial codes and indexing eye data
    times_ms = cellfun(@(q) round(q(1)*1000:q(2)*1000), num2cell(tbl1.time_sec,2), 'uni', 0);
    trialStarts = cellfun(@(q,r) find(q == round(r(r(:,2)==1,3)*1000)), times_ms, tbl1.trialcodes, 'uni', 0);
    eventCodes = tbl1.trialcodes; eventCodes = vertcat(eventCodes{:});
    eventCodes = num2cell(sort(unique(eventCodes(:,2))))';

    eventNames = convertBetween_eventCodes_eventNames(eventCodes);
    trialMarkers = cellfun(@(t) cellfun(@(q,r) find(ismember(q,round(r(r(:,2)==t,3)*1000))), times_ms, tbl1.trialcodes, 'uni', 0), eventCodes, 'uni', 0)';
    trialMarkers = horzcat(trialMarkers{:});
    trialMarkers(cellfun('isempty',trialMarkers)) = {NaN};

    for m=1:length(eventNames)
        if sum(cellfun(@(q) size(q,1)>1, trialMarkers(:,m), 'uni', 1))>0 || sum(cellfun(@(q) size(q,2)>1, trialMarkers(:,m), 'uni', 1))>0
            tbl.(eventNames{m}) = trialMarkers(:,m);
        else
            tbl.(eventNames{m}) = cell2mat(trialMarkers(:,m));
        end
    end

    tbl.block = tbl1.block; 
    tbl.params = tbl1.params;

    eyePos = cellfun(@(x,y) filterEyeTraces_EyeLink(x(:,y:end),'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',84,'PLOT_TRIAL',false), tbl1.eyes, trialStarts, 'uni', 0);
    eyeVel = cellfun(@(q) calcDerivative_eyeTraces(q), cellfun(@(x,y) filterEyeTraces_EyeLink(x(:,y:end),'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',40,'PLOT_TRIAL',false), tbl1.eyes, trialStarts, 'uni', 0), 'uni', 0);
    eyeAcc = cellfun(@(q) calcDerivative_eyeTraces(q), eyeVel, 'uni', 0);

    tbl.eyePos = eyePos; tbl.eyeVel = eyeVel; tbl.eyeAcc = eyeAcc;
    tbl.pupil = tbl1.pupil; 
    tbl.spiketimes = tbl1.spiketimes; tbl.net_labels = tbl1.net_labels;

    %%%%%%%%%%%%% task-specific re-arranging and calculations %%%%%%%%%%%%%
    if contains(task_name, 'mdir')
        tbl(:, 'MEM_GUIDED_SACC') = [];
        tbl.distance = cellfun(@(q) round(pix2deg(q,tbl(1,:).params.block.screenDistance,tbl(1,:).params.block.pixPerCM)), num2cell(tbl.distance), 'uni', 1);

        if all(ismember({'targetOnsetDelay', 'delay'}, tbl.Properties.VariableNames))
            tbl.fixDuration = tbl.targetOnsetDelay+tbl(1,:).params.block.targetDuration+tbl.delay;

            tbl = movevars(tbl,{'targetOnsetDelay','delay','fixDuration'},'Before','result');
        end
    elseif contains(task_name, 'purs')
        % Define the columns to replace and their new names
        cols_to_replace = {'TARG_ON', 'BROKE_TARG', 'TARG_OFF'};
        new_names = {'PURSUIT_TARG', 'BROKE_PURSUIT', 'PURSUIT_TARG_OFF'};
        
        % Loop through each column to check and replace
        for i = 1:numel(cols_to_replace)
            if ismember(cols_to_replace{i}, tbl.Properties.VariableNames)
                tbl.Properties.VariableNames{ismember(tbl.Properties.VariableNames, cols_to_replace{i})} = new_names{i};
            end
        end

        [pursuitOnsets, rxnTimes, msFlags, csFlags] = deal(zeros(height(tbl), 1));
        for t = 1:height(tbl)
            if isequal(tbl.result(t),"CORRECT")
                [~,radPos] = cart2pol(tbl.eyePos{t}(1,:),tbl.eyePos{t}(2,:));
                [~,radVel] = cart2pol(tbl.eyeVel{t}(1,:),tbl.eyeVel{t}(2,:));
                [pursuitOnset,rxnTime,msFlag,csFlag] = detect_pursuitOnset(radPos,radVel,tbl.PURSUIT_TARG(t),tbl(t,:).params.block.reactionTime,tbl.pursuitSpeed(t),tbl.angle(t),tbl.jump(t),0);
            else
                pursuitOnset = NaN;
                rxnTime = NaN; 
                msFlag = NaN;
                csFlag = NaN;
            end
            pursuitOnsets(t) = pursuitOnset;
            rxnTimes(t) = rxnTime;
            msFlags(t) = msFlag;
            csFlags(t) = csFlag;
        end

        tbl.pursuitOnset = pursuitOnsets;
        tbl.rxnTime = rxnTimes;
        tbl.msFlag = msFlags;
        tbl.csFlag = csFlags;

        tbl = movevars(tbl,{'pursuitOnset','rxnTime','msFlag','csFlag'},'Before','result');

    elseif contains(task_name, 'rfmp')
        tbl = tbl(cellfun(@(q) sum(isnan(q))==0, tbl.STIM_ON, 'uni', 1),:);
        tbl.STIM_ON(tbl.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl.STIM_ON(tbl.result~='CORRECT'), 'uni', 0);
        tbl.conditions(tbl.result~="CORRECT") = cellfun(@(q) q(1:end-1), tbl.conditions(tbl.result~='CORRECT'), 'uni', 0);
    end

    if ismember('emptyCnd', tbl.Properties.VariableNames)
        tbl.emptyCnd = [];
    end

end
