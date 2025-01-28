function dat_all = format_datTrials(nev1,out,eye_channel_labels)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    Fs = double(out.hdr.Fs);
    nEpochs = size(out.hdr.timeStamps,2);
    
    starttrial = 1;
    endtrial = 255;
    
    % {Eye_HE, Eye_VE, DIODE, PUPIL}
    eye_channels = find(ismember(out.hdr.label, eye_channel_labels));
    
    % determine if nev is an array of struct, if struct pull out nev
    if isequal(class(nev1),'struct')
        nev = [nev1.nev nev1.net_labels'];
        %chan_units = unique(nev(nev(:,1) ~= 0 & nev(:,2) ~=0, 1:2), 'rows');
        spike_sort = true;
    else
        nev = nev1;
        spike_sort = false;
    end

    spk_channels = sort(unique(nev(nev(:,1) ~= 0, 1)));
    
    dat_all = [];
    past_epochEnd = 0;
    block = 1;
    for epoch=1:nEpochs
        % pull out data for this epoch
        epochStart = out.hdr.timeStamps(1,epoch); % samp
        epochEnd = out.hdr.timeStamps(2,epoch); % samp
    
        nsStartTime = double(epochStart / Fs); % sec
        nsEndTime = double(epochEnd / Fs) + 0.3; % sec
    
        epochDiff = epochEnd - epochStart;
        epochStart_samp = past_epochEnd + 1;
        epochEnd_samp = (epochStart_samp + epochDiff)-1;
    
        this_nev = nev(nev(:,3)>=nsStartTime & nev(:,3)<=nsEndTime,:);
        this_ns5 = out.data(:,epochStart_samp:epochEnd_samp);

        diginnevind = find(this_nev(:,1)==0);
        digcodes = this_nev(diginnevind,:);
    
        trialstartindstemp = (find(digcodes(:,2)==starttrial));
        trialstartinds = diginnevind(trialstartindstemp);
        trialstarts = this_nev(trialstartinds,3);
    
        trialendindstemp = (find(digcodes(:,2)==endtrial));
        trialendinds = diginnevind(trialendindstemp);
        trialends = this_nev(trialendinds,3);
    
        [trialstarts, trialends,trialstartgood,trialendgood] = detectMissingStartEndCode(trialstarts,trialends);
        trialstartinds = trialstartinds(trialstartgood);
        trialendinds = trialendinds(trialendgood);
    
        if length(trialstarts)~=length(trialends) || sum((trialends-trialstarts)<0)
            % fix it
            if sum(trialstarts(1:end-1)>=trialends)==0
                trialstarts = trialstarts(1:end-1);
            end
        end
    
        trialstarts_samp = round(trialstarts*Fs) - epochStart;
        trialends_samp = round(trialends*Fs) - epochStart;
        past_epochEnd = epochEnd_samp;
    
        % Get session initial params
        
        predatcodes = digcodes(digcodes(:,3)<trialstarts(1),:);
        tempdata.text = char(predatcodes(predatcodes(:,2)>=256 & predatcodes(:,2)<512,2)-256)';
        if ~isempty(tempdata.text)
    
            tempdata = getDatParams(tempdata);
    
            dat = struct();
    
            %% Make Struct
            for n = 1:length(trialstarts)
                if mod(n,100) == 0
                    fprintf('Processed nev for %i trials of %i...\n',n,length(trialstarts));
                end
                dat(n).block = block;
                dat(n).time_sec = [trialstarts(n) trialends(n)];
                this_trial = this_nev(trialstartinds(n):trialendinds(n),:);
                trialdig = this_trial(this_trial(:,1)==0,:);
                dat(n).text = char(trialdig(trialdig(:,2)>=256 & trialdig(:,2)<512,2)-256)';
                dat(n).trialcodes = trialdig(trialdig(:,2)<256 | (trialdig(:,2)>=1000 & trialdig(:,2)<=32000),:);
    
                event = uint32(trialdig);
                dat(n).result = event(event(:,2)>=160 & event(:,2)<=165,2);
                if isempty(dat(n).result)
                    dat(n).result = event(event(:,2)>=150 & event(:,2)<=158,2);
                end
                if(isempty(dat(n).result))
                    dat(n).result = NaN;
                end
    
                dat(n).params.block = tempdata.params.trial;
                if n<length(trialstarts) && trialstartinds(n+1)- trialendinds(n)>1
                    bt = this_nev(trialendinds(n)+1:trialstartinds(n+1)-1,:);
                    btdig = bt(bt(:,1)==0,:);
                    if sum(find(btdig(:,2)>=256 & btdig(:,2)<512))> 0
                        tempdata.text = char(btdig(btdig(:,2)>=256 & btdig(:,2)<512,2)-256)';
                        tempdata = getDatParams(tempdata);
                        block = block + 1;
                    end
                end
    
                % eye data
                eyes = this_ns5(eye_channels([1,2,4]),trialstarts_samp(n):trialends_samp(n));
                eyes_1khz = downsample(eyes',30)';
                [eyedeg, ~] = eye2deg(eyes_1khz(1:2,:), dat(n).params);
    
                pupil_1kHz = eyes_1khz(3,:);
    
                dat(n).eyes = eyedeg;
                dat(n).pupil = pupil_1kHz;
                dat(n).diode = this_ns5(eye_channels(3),trialstarts_samp(n):trialends_samp(n));
                dat(n).ns5_30kHz = this_ns5(:,trialstarts_samp(n):trialends_samp(n));
    
                % spikes
                spks = this_trial(ismember(this_trial(:,1), spk_channels),:);
                spks_byChan = cell(1,size(spk_channels,1));
                if spike_sort
                    % [netLabels_byChan,waveforms_byChan] = deal(cell(1,size(spk_channels,1)));
                    netLabels_byChan = cell(1,size(spk_channels,1));
                end
                for u = 1:length(spks_byChan)
                    spks_byChan{u} = (spks(ismember(spks(:, 1), spk_channels(u)), 3)') - trialstarts(n);
                    if spike_sort
                        netLabels_byChan{u} = (spks(ismember(spks(:, 1), spk_channels(u)), 4)');
                        % waveforms_byChan{u} = (spks(ismember(spks(:, 1), spk_channels(u)), 5:end)');
                    end
                end
                dat(n).spiketimes = spks_byChan;
                if spike_sort
                    dat(n).net_labels = netLabels_byChan;
                    % dat(n).waveforms = waveforms_byChan;
                end
    
            end
            dat = getDatParams(dat);
        end
    
        % Concatenate the new structure to the array
        dat_all = [dat_all, dat];
    
    end
end