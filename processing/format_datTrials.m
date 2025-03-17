function [dat_all,epochEnd] = format_datTrials(nev1, out_ns5, eye_channel_labels, neural_channels)
    % format_datTrials - Processes neural and behavioral data for multiple trials, 
    % extracts eye and spike data, and formats the information into a structured array.
    %
    % This function processes and organizes trial-based data from neural event (nev) and neural signal (ns5) files. 
    % It extracts relevant trial-specific information such as eye position, velocity, and pupil data, as well as spike 
    % times for specified neural channels. The function is capable of handling various data epochs and organizes the data 
    % into a structured format that includes trial codes, trial results, parameters, and neural spikes.
    %
    % The function also supports neural spike sorting and handles missing or corrupted trial start/end codes. 
    % The data is downsampled and formatted for further analysis.
    %
    %%%% Required inputs: %%%
    %   nev1               -   nev file containing event data (e.g., neural threshold crossings and digital codes)
    %   out_ns5            -   ns5 file containing the raw 30kHz data (e.g., eye data, pupil, diode, raw neural signals)
    %   eye_channel_labels -   Cell array of labels for eye movement channels (Eye_HE', 'Eye_VE', 'DIODE', 'PUPIL')
    %                          e.g. {'10241', '10242', '10243', '10244'}
    %   neural_channels    -   Array of channel IDs to extract neural spike data from (e.g., [1,2,3,4,5])
    %
    %%%% Outputs: %%%
    %   dat_all         -    Array of structs with the following fields:
    %       - block: The block number of the trial
    %       - time_sec: Start and end times of the trial in seconds
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
    %%%% Example usage: %%%
    %   dat_all = format_datTrials(nev1, out_ns5, {'10241', '10242', '10243', '10244'}, [0,1,2,3,4,5,6,7,8,9])
    %
    % This will process the nev1 and out_ns5 files, extract eye and neural data, 
    % and return the data in the structured array `dat_all`, with separate trial information for each epoch.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fs = double(out_ns5.hdr.Fs); % Sampling frequency
    nEpochs = size(out_ns5.hdr.timeStamps, 2); % Number of epochs
    
    starttrial = 1; % Trial start code
    endtrial = 255; % Trial end code
    
    % Identify eye movement channels (Eye_HE, Eye_VE, DIODE, PUPIL)
    if ~isempty(eye_channel_labels)
        eye_channels = find(ismember(out_ns5.hdr.label, eye_channel_labels));
    end
    
    % Determine if nev is an array of struct; if so, extract nev data
    if isequal(class(nev1), 'struct')
        nev = [nev1.nev nev1.net_labels']; % Combine event and neural network labels
        spike_sort = true; % Flag for spike sorting
    else
        nev = nev1;
        spike_sort = false; % No spike sorting if nev is not a struct
    end
    
    dat_all = []; % Initialize output structure array
    past_epochEnd = 0; % Keep track of the end of the previous epoch
    block = 1; % Block number for trial grouping
    
    % Loop through each epoch
    for epoch = 1:nEpochs
        % Extract data for the current epoch
        epochStart = out_ns5.hdr.timeStamps(1, epoch); % Epoch start time (samples)
        epochEnd = out_ns5.hdr.timeStamps(2, epoch); % Epoch end time (samples)
    
        nsStartTime = double(epochStart / Fs); % Convert to seconds
        nsEndTime = double(epochEnd / Fs) + 0.3; % Add a small buffer for the end time
    
        epochDiff = epochEnd - epochStart;
        epochStart_samp = past_epochEnd + 1;
        epochEnd_samp = (epochStart_samp + epochDiff) - 1;
    
        % Extract trial data within the current epoch time range
        this_nev = nev(nev(:, 3) >= nsStartTime & nev(:, 3) <= nsEndTime, :);
        ns5_rng = epochStart_samp:epochEnd_samp; % Range of samples in the epoch

        % Extract digital event codes
        diginnevind = find(this_nev(:, 1) == 0);
        digcodes = this_nev(diginnevind, :);
    
        % Find trial start and end indices
        trialstartindstemp = find(digcodes(:, 2) == starttrial);
        trialstartinds = diginnevind(trialstartindstemp);
        trialstarts = this_nev(trialstartinds, 3); % Trial start times

        trialendindstemp = find(digcodes(:, 2) == endtrial);
        trialendinds = diginnevind(trialendindstemp);
        trialends = this_nev(trialendinds, 3); % Trial end times
    
        % Detect missing start/end codes and handle missing data
        [trialstarts, trialends, trialstartgood, trialendgood] = detectMissingStartEndCode(trialstarts, trialends);
        trialstartinds = trialstartinds(trialstartgood);
        trialendinds = trialendinds(trialendgood);
    
        % Handle mismatched start/end times
        if length(trialstarts) ~= length(trialends) || sum((trialends - trialstarts) < 0)
            if sum(trialstarts(1:end-1) >= trialends) == 0
                trialstarts = trialstarts(1:end-1); % Trim the last trial start
            end
        end
    
        % Convert trial start and end times to samples
        trialstarts_samp = round(trialstarts * Fs) - epochStart;
        trialends_samp = round(trialends * Fs) - epochStart;
        past_epochEnd = epochEnd_samp;
    
        % Get session initial parameters
        predatcodes = digcodes(digcodes(:, 3) < trialstarts(1), :);
        tempdata.text = char(predatcodes(predatcodes(:, 2) >= 256 & predatcodes(:, 2) < 512, 2) - 256)';
        if ~isempty(tempdata.text)
    
            tempdata = getDatParams(tempdata); % Get parameters associated with the data
            if ~isempty(eye_channel_labels)
                dat = repmat(struct(...
                    'block', [], ...
                    'time_sec', [], ...
                    'text', '', ...
                    'trialcodes', [], ...
                    'result', NaN, ...
                    'params', struct(), ...
                    'eyes', [], ...
                    'pupil', [], ...
                    'diode', []), length(trialstarts), 1);
            else
                dat = repmat(struct(...
                    'block', [], ...
                    'time_sec', [], ...
                    'text', '', ...
                    'trialcodes', [], ...
                    'result', NaN, ...
                    'params', struct()), length(trialstarts), 1);
            end

            %% Loop through trials and organize data
            for n = 1:length(trialstarts)
                if mod(n, 100) == 0
                    fprintf('Processed nev for %i trials of %i...\n', n, length(trialstarts));
                end
                dat(n).block = block;
                dat(n).time_sec = [trialstarts(n) trialends(n)]; % Store trial start and end times
                this_trial = this_nev(trialstartinds(n):trialendinds(n), :);
                trialdig = this_trial(this_trial(:, 1) == 0, :);
                dat(n).text = char(trialdig(trialdig(:, 2) >= 256 & trialdig(:, 2) < 512, 2) - 256)';
                dat(n).trialcodes = trialdig(trialdig(:, 2) < 256 | (trialdig(:, 2) >= 1000 & trialdig(:, 2) <= 32000), :);
    
                % Extract trial result
                event = uint32(trialdig);
                dat(n).result = event(event(:, 2) >= 160 & event(:, 2) <= 165, 2);
                if isempty(dat(n).result)
                    dat(n).result = event(event(:, 2) >= 150 & event(:, 2) <= 158, 2);
                end
                if isempty(dat(n).result)
                    dat(n).result = NaN;
                end
    
                % Extract block parameters
                blockParams = tempdata.params.trial;
                if isfield(blockParams, 'reactionTime')
                    blockParams.crossingTime = blockParams.reactionTime; % Copy reaction time to crossing time
                    blockParams = rmfield(blockParams, 'reactionTime'); % Remove old field
                end
                                
                dat(n).params.block = blockParams;
                
                % Check if next trial is in a new block
                if n < length(trialstarts) && trialstartinds(n+1) - trialendinds(n) > 1
                    bt = this_nev(trialendinds(n) + 1:trialstartinds(n + 1) - 1, :);
                    btdig = bt(bt(:, 1) == 0, :);
                    if sum(find(btdig(:, 2) >= 256 & btdig(:, 2) < 512)) > 0
                        tempdata.text = char(btdig(btdig(:, 2) >= 256 & btdig(:, 2) < 512, 2) - 256)';
                        tempdata = getDatParams(tempdata);
                        block = block + 1; % Move to the next block
                    end
                end

                dat(n).ns5_samps = ns5_rng([trialstarts_samp(n),trialends_samp(n)]);

                % Extract and process eye data
                if ~isempty(eye_channel_labels)
                    eyes = out_ns5.data(eye_channels([1, 2, 3, 4]), ns5_rng(trialstarts_samp(n):trialends_samp(n)));
                    eyes_1khz = downsample(eyes', 30)'; % Downsample to 1 kHz
                    [eyedeg, ~] = eye2deg(eyes_1khz(1:2, :), dat(n).params); % Convert to degrees
        
                    % Store eye, pupil, and diode data
                    dat(n).eyes = eyedeg;
                    dat(n).pupil = eyes_1khz(4, :);
                    dat(n).diode = eyes_1khz(3, :);
                end
                
                % Process neural spikes, if applicable
                if ~isempty(neural_channels)
                    spks = this_trial(ismember(this_trial(:, 1), neural_channels), :);
                    spks_byChan = cell(1, size(neural_channels, 1));
                    if spike_sort
                        [netLabels_byChan, waveforms_byChan] = deal(cell(1, size(neural_channels, 1)));
                    end
                    for u = 1:length(spks_byChan)
                        spks_byChan{u} = ((spks(ismember(spks(:, 1), neural_channels(u)), 3)') - trialstarts(n)) .* 1000;
                        if spike_sort
                            netLabels_byChan{u} = (spks(ismember(spks(:, 1), neural_channels(u)), 4)');
                            waveforms_byChan{u} = (spks(ismember(spks(:, 1), neural_channels(u)), 5:end)');
                        end
                    end
                    dat(n).spiketimes = spks_byChan; % Store spike times
                    if spike_sort
                        dat(n).net_labels = netLabels_byChan; % Store spike sorting labels
                    end
                end
    
            end
            dat = getDatParams(dat); % Final parameter extraction
        end
    
        % Concatenate the new structured data to the main array
        dat_all = [dat_all; dat];
    end
end