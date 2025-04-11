function eventsOut = convertBetween_eventCodes_eventNames(eventsIn)
    % convertBetween_eventCodes_eventNames - Converts event codes to event names or vice versa.
    %
    % This function takes in an array of event codes (numeric values) or event names (strings/categoricals), 
    % and returns the corresponding array of event names or event codes.
    % It uses a predefined mapping of event codes to event names and vice versa, allowing for 
    % easy conversion between these representations of events.
    %
    %%%% Inputs: %%%
    %   eventsIn    -   An array of event codes (numeric or cell array of strings/categoricals).
    %
    %%%% Outputs: %%%
    %   eventsOut   -   An array of event names (if input is event codes) or event codes (if input is event names).
    %                   The output format matches the input format (numeric or string).
    %
    % The mapping between event codes and names is defined in the 'codes' structure.
    %
    %%%% Example usage: %%%
    %   eventsOut = convertBetween_eventCodes_eventNames([1, 2, 3])
    %   % Returns event names corresponding to codes 1, 2, and 3.
    %
    %   eventsOut = convertBetween_eventCodes_eventNames({'START_TRIAL', 'FIX_ON'})
    %   % Returns event codes corresponding to 'START_TRIAL' and 'FIX_ON'.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define the mapping between event codes and event names
    codes = struct();

    % Trial start and end codes
    codes.ALIGN_PULSE = 0;
    codes.START_TRIAL = 1;
    codes.BCI_END = 249;
    codes.BACKGROUND_PROCESS_TRIAL = 250;
    codes.SHOWEX_TIMINGERROR = 251;
    codes.BCI_ABORT = 252;
    codes.ALIGN = 253;
    codes.SHOWEX_ABORT = 254;
    codes.END_TRIAL = 255;

    % Stimulus / trial event codes
    codes.FIX_ON = 2;
    codes.FIX_OFF = 3;
    codes.FIX_MOVE = 4;
    codes.REWARD = 5;
    codes.DIODE_ON = 6;
    codes.DIODE_OFF = 7;
    codes.NO_STIM = 9;
    codes.STIM_ON = 10;
    codes.STIM1_ON = 11;
    codes.STIM2_ON = 12;
    codes.STIM3_ON = 13;
    codes.STIM4_ON = 14;
    codes.STIM5_ON = 15;
    codes.STIM6_ON = 16;
    codes.STIM7_ON = 17;
    codes.STIM8_ON = 18;
    codes.STIM9_ON = 19;
    codes.STIM10_ON = 20;
    codes.STIM_OFF = 40;
    codes.STIM1_OFF = 41;
    codes.STIM2_OFF = 42;
    codes.STIM3_OFF = 43;
    codes.STIM4_OFF = 44;
    codes.STIM5_OFF = 45;
    codes.STIM6_OFF = 46;
    codes.STIM7_OFF = 47;
    codes.STIM8_OFF = 48;
    codes.STIM9_OFF = 49;
    codes.STIM10_OFF = 50;
    codes.TARG_ON = 70;
    codes.TARG1_ON = 71;
    codes.TARG2_ON = 72;
    codes.TARG3_ON = 73;
    codes.TARG4_ON = 74;
    codes.TARG5_ON = 75;
    codes.TARG6_ON = 76;
    codes.TARG7_ON = 77;
    codes.TARG8_ON = 78;
    codes.TARG9_ON = 79;
    codes.TARG10_ON = 80;
    codes.TARG_OFF = 100;
    codes.TARG1_OFF = 101;
    codes.TARG2_OFF = 102;
    codes.TARG3_OFF = 103;
    codes.TARG4_OFF = 104;
    codes.TARG5_OFF = 105;
    codes.TARG6_OFF = 106;
    codes.TARG7_OFF = 107;
    codes.TARG8_OFF = 108;
    codes.TARG9_OFF = 109;
    codes.TARG10_OFF = 110;
    codes.CHOICE0 = 120;
    codes.CHOICE1 = 121;
    codes.CHOICE2 = 122;
    codes.CHOICE3 = 123;
    codes.CHOICE4 = 124;
    codes.CHOICE5 = 125;
    codes.CHOICE6 = 126;
    codes.CHOICE7 = 127;
    codes.CHOICE8 = 128;
    codes.CHOICE9 = 129;

    % UStim codes
    codes.USTIM_ON = 130;
    codes.USTIM_OFF = 131;

    % Sound codes
    codes.SOUND_ON = 132;
    codes.SOUND_OFF = 133;
    codes.SOUND_CHANGE = 134;
    codes.CURSOR_ON = 135;
    codes.CURSOR_OFF = 136;

    % Behavior codes
    codes.FIXATE = 140; % attained fixation
    codes.SACCADE = 141; % initiated saccade
    codes.CURSOR_POS = 142; % indicates next codes will define cursor position
    codes.BCI_CURSOR_POS = 143; % indicates next codes will define cursor position from BCI

    % Trial outcome codes
    codes.CORRECT = 150; % Independent of whether reward is given
    codes.IGNORED = 151; % Never fixated or started trial
    codes.BROKE_FIX = 152; % Left fixation before trial complete
    codes.WRONG_TARG = 153; % Chose wrong target
    codes.BROKE_TARG = 154; % Left target fixation before required time
    codes.MISSED = 155; % for a detection task
    codes.FALSEALARM = 156;
    codes.NO_CHOICE = 157; % saccade to non-target / failure to leave fix window
    codes.WITHHOLD = 158; % correctly-withheld response
    codes.ACQUIRE_TARG = 159; % Acquired the target
    codes.FALSE_START = 160; % left too early
    codes.BCI_CORRECT = 161; % BCI task performed correctly
    codes.BCI_MISSED = 162; % BCI task performed incorrectly
    codes.CORRECT_REJECT = 163;
    codes.LATE_CHOICE = 164;
    codes.BROKE_TASK = 165;
    codes.PURSUIT_TARG = 166;
    codes.BROKE_PURSUIT = 167;
    codes.PURSUIT_TARG_ON = 31791;
    codes.PURSUIT_TARG_OFF = 12697;

    % Block event codes
    codes.SACCADE_BLOCK = 1000;
    codes.PURSUIT_BLOCK = 1001;

    % Guided saccade event codes
    codes.VIS_GUIDED_SACC = 2001;
    codes.MEM_GUIDED_SACC = 2002;
    codes.DVIS_GUIDED_SACC = 2003;

    % NaN code
    codes.NaN = 666;

    % Convert the structure to cell arrays for easy access
    allCodes = struct2cell(codes);
    allEvents = fieldnames(codes);

    %% CONVERSION LOGIC
    % Check if input events are numeric codes or event names and perform conversion
    if isequal(class(eventsIn), 'double') || isequal(class(eventsIn), 'uint32')
        eventsIn = num2cell(eventsIn);  % Ensure eventsIn is a cell array
    end

    % Case 1: If input is event names (strings), convert them to event codes
    if isequal(class(eventsIn{1}), 'char')
        [~, index] = ismember(eventsIn, allEvents);
        eventsOut = allCodes(index);
        
    % Case 2: If input is event names (string/categorical), convert to event codes
    elseif isequal(class(eventsIn{1}), 'string') || isequal(class(eventsIn{1}), 'categorical')
        eventsIn = cellfun(@(q) char(q), eventsIn, 'uni', 0); % Convert string/categorical to char
        [~, index] = ismember(eventsIn, allEvents);
        eventsOut = allCodes(index);
        
    % Case 3: If input is numeric event codes, convert them to event names
    elseif isequal(class(eventsIn{1}), 'double') || isequal(class(eventsIn{1}), 'uint32')
        eventsIn = cellfun(@(q) cast(q, 'double'), eventsIn, 'uni', 0); % Ensure double type
        eventsIn = cellfun(@(q) q(1), eventsIn, 'uni', 0); % Handle case where two results are reported
        eventsIn = cell2mat(eventsIn);
        eventsIn(isnan(eventsIn)) = 666; % Handle NaN as 'NaN' event code
        [~, index] = ismember(eventsIn, cell2mat(allCodes));
        eventsOut = allEvents(index); % Output event names
    end
end