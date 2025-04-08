function raw_all = extract_rawData(nev1,out,channels)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    Fs = double(out.hdr.Fs);
    nEpochs = size(out.hdr.timeStamps,2);
    
    starttrial = 1;
    endtrial = 255;

    neural_inds = find(ismember(out.hdr.label,string(1:500)));
    
    % determine if nev is an array of struct, if struct pull out nev
    if isequal(class(nev1),'struct')
        nev = [nev1.nev nev1.net_labels'];
    else
        nev = nev1;
    end
    
    raw_all = []; %waves_all = [];
    past_epochEnd = 0;
    for epoch=1:nEpochs
        % pull out data for this epoch
        epochStart = (out.hdr.timeStamps(1,epoch))./(30000/Fs); % samp
        epochEnd = (out.hdr.timeStamps(2,epoch))./(30000/Fs); % samp
    
        nsStartTime = double(epochStart / Fs); % sec
        nsEndTime = double(epochEnd / Fs) + (Fs/100000); % sec
    
        epochDiff = epochEnd - epochStart;
        epochStart_samp = past_epochEnd + 1;
        epochEnd_samp = (epochStart_samp + epochDiff)-1;
    
        this_nev = nev(nev(:,3)>=nsStartTime & nev(:,3)<=nsEndTime,:);
        ns5_rng = epochStart_samp:epochEnd_samp;

        diginnevind = find(this_nev(:,1)==0);
        digcodes = this_nev(diginnevind,:);
    
        trialstartindstemp = (find(digcodes(:,2)==starttrial));
        trialstartinds = diginnevind(trialstartindstemp);
        trialstarts = this_nev(trialstartinds,3);
    
        trialendindstemp = (find(digcodes(:,2)==endtrial));
        trialendinds = diginnevind(trialendindstemp);
        trialends = this_nev(trialendinds,3);
    
        [trialstarts, trialends,~,~] = detectMissingStartEndCode(trialstarts,trialends);
    
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
            rawData = cell(length(trialstarts),1);
    
            %% Make Struct
            for n = 1:length(trialstarts)
                if mod(n,100) == 0
                    fprintf('Processed nev for %i trials of %i...\n',n,length(trialstarts));
                end

                neural_data = out.data(neural_inds,ns5_rng(trialstarts_samp(n):trialends_samp(n)));
                rawData{n} = neural_data(channels,:);

            end
        end
    
        % Concatenate the new structure to the array
        raw_all = [raw_all; rawData];

    end
end