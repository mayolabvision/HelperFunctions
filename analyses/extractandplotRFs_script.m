%% LOAD THE DATA FILE
load kendra_scrappy_0132a_g0_unleashed.mat


%% SCRIPT TO EXTRACT RF INFO FROM DATA TABLE (USING RELATIVE SPIKE TIMES STORED IN 'spiketimes_imec1' in S.rfmp1.data)

%GET STIMULUS info
T = S.rfMapping_dots_0001.tbl;
all_stimuli = []
%get all unique conditions
for rr = 1:size(T,1) %loop over rows of the table
    if isnan(T{rr,'BROKE_FIX'}) %if no broken fixation
        stims = T{rr,'conditions'}{1};%pull stimuli
        for ss = 1:size(stims,1) %loop over stimuli
            all_stimuli = [all_stimuli; stims{ss,1}]; %concatenate stimuli
        end
    end
end

unique_xvals = unique(all_stimuli(:,1));
unique_yvals = unique(all_stimuli(:,1));
[Xgrid,Ygrid] = meshgrid(unique_xvals, flipdim(unique_yvals,1)); %grids that represent X and Y coordinates on screen

% collect addresses of where each stimulus occurs (registered to the Xgrid,Ygrid matrices)
Tstimloc = cell(size(Xgrid)); %make empty cell array to store matrix of data table locations for these corresponding stimuli [[trial epochidx]; ...]

for rr = 1:size(T,1) %loop over trials
    if isnan(T{rr,'BROKE_FIX'}) %if no broken fixation
        stims = T{rr,'conditions'}{1};%pull stimuli
        for ss = 1:size(stims,1) %loop over stimuli
            thisstim = stims{ss,1}; %pull this stimulus
            [idxRow,idxCol] = find((Xgrid==thisstim(1)) & (Ygrid==thisstim(2))); %find where this stimulus sits in the stim grid matrices
            
            locations_so_far = Tstimloc{idxRow, idxCol}; %pull the stimulus locations for this stimulus we've stored up to now
            updated_locations = [locations_so_far; [rr ss]]; %append the current ROW and EPOCH to the location list
            Tstimloc{idxRow, idxCol} = updated_locations;
        end
    end
    
end

%set time bin offsets (from STIM_ON time)
twin=50; %time bin window (in ms) to count spikes over
offsets(:,1) = [0:10:250]';
offsets(:,2) = offsets(:,1)+twin;




% Build up RF maps and store for each probe
IMEC1 = [];

imec_str = 'spiketimes_imec1'; %define the probe data to pull spikes from
n_units = size(T{1,imec_str}{1},2); %get number N of units

for uu = 1:n_units %loop over units
    disp(['mapping unit ' num2str(uu)])
    %tic
    avg_rf_map = nan(size(Xgrid,1), size(Xgrid,2), size(offsets,1)); %will store average spike counts [Xpos x Ypos X Timebin]
    rf_map_cell = cell(size(Xgrid,1), size(Xgrid,2), size(offsets,1)); %will store spike counts for each presentation in a cell array
    
    for rr=1:size(Tstimloc,1) %rows in the location matrix
        for cc=1:size(Tstimloc,2) %columns in the location matrix
            %disp(['XY = [' num2str(Xgrid(rr,cc)) ',' num2str(Ygrid(rr,cc)) ']'])
            
            locations = Tstimloc{rr,cc};
            for ll = 1:size(locations,1) %loop over locations (aka data address of each stim presentation)
                thisloc = locations(ll,:); %location is [trialnum epochnum]
                stim_on = T{thisloc(1),'STIM_ON'}{1}(thisloc(2)); %relative time of stimulus on in the trial (in milliseconds)
                trialspikes = round(T{thisloc(1),imec_str}{1}{1,uu}); %relative spike times for this neuron on this trial (in milliseconds)
                
                for oo = 1:size(offsets,1) %loop over time bins to consider
                   spike_rng = stim_on + offsets(oo,:); % sets relative time window to count spikes in
                   spike_ct = sum(trialspikes>=spike_rng(1) & trialspikes<spike_rng(2));
                   
                   previous_counts = rf_map_cell{rr,cc,oo};
                   new_counts = [previous_counts; spike_ct];
                   rf_map_cell{rr,cc,oo} = new_counts;
                   
                end
                
            end
        end
    end
    
    avg_rf_map = cellfun(@nanmean, rf_map_cell);
    avg_rf_map = avg_rf_map*(1000/twin); %convert avg rf map to Hz firing rates
    
    IMEC1.unit(uu).avg_rf_map = avg_rf_map;
    IMEC1.unit(uu).rf_map_cell = rf_map_cell;   
    %toc
    
end










%% Code as above, but using S.Kilosort timestamps for the spike times, and registering with S.rfmmp1.data{:,'ns5_samp'} times (presumably trial start/end in ns5 timestamp format) and S.rfmp1.data{:,'STIM_ON'} times (presumably relative time within the trial for each stimulus onset in milliseconds)

%GET STIMULUS info
T = S.rfmp1.data;
all_stimuli = []
%get all unique conditions
for rr = 1:size(T,1) %loop over rows of the table
    if isnan(T{rr,'BROKE_FIX'}) %if no broken fixation
        stims = T{rr,'conditions'}{1};%pull stimuli
        for ss = 1:size(stims,1) %loop over stimuli
            all_stimuli = [all_stimuli; stims{ss,1}]; %concatenate stimuli
        end
    end
end

unique_xvals = unique(all_stimuli(:,1));
unique_yvals = unique(all_stimuli(:,1));
[Xgrid,Ygrid] = meshgrid(unique_xvals, flipdim(unique_yvals,1)); %grids that represent X and Y coordinates on screen

% collect addresses of where each stimulus occurs (registered to the Xgrid,Ygrid matrices)
Tstimloc = cell(size(Xgrid)); %make empty cell array to store matrix of data table locations for these corresponding stimuli [[trial epochidx]; ...]

for rr = 1:size(T,1) %loop over trials
    if isnan(T{rr,'BROKE_FIX'}) %if no broken fixation
        stims = T{rr,'conditions'}{1};%pull stimuli
        for ss = 1:size(stims,1) %loop over stimuli
            thisstim = stims{ss,1}; %pull this stimulus
            [idxRow,idxCol] = find((Xgrid==thisstim(1)) & (Ygrid==thisstim(2))); %find where this stimulus sits in the stim grid matrices
            
            locations_so_far = Tstimloc{idxRow, idxCol}; %pull the stimulus locations for this stimulus we've stored up to now
            updated_locations = [locations_so_far; [rr ss]]; %append the current ROW and EPOCH to the location list
            Tstimloc{idxRow, idxCol} = updated_locations;
        end
    end
    %disppercent(rr/size(T,1))
end

%set time bin offsets (from STIM_ON time)
twin=50; %time bin window (in ms) to count spikes over
offsets(:,1) = [0:10:250]';
offsets(:,2) = offsets(:,1)+twin;

% Build up RF maps and store for each probe IMEC 0
imec_num=0;
eval(['IMEC' num2str(imec_num) ' = [];']);


imec_str = ['spiketimes_imec' num2str(imec_num)]; %define the probe data to pull spikes from
n_units = size(T{1,imec_str}{1},2); %get number N of units
clusters = S.kilosort(imec_num+1).clusters{:,'cluster'}; %get unique cluster #s
spike_clusters = S.kilosort(imec_num+1).spike_clusters;%get cluster #s for every spike

for uu = 1:n_units %loop over units
    disp(['mapping unit ' num2str(uu)])
    tic
    avg_rf_map = nan(size(Xgrid,1), size(Xgrid,2), size(offsets,1)); %will store average spike counts [Xpos x Ypos X Timebin]
    rf_map_cell = cell(size(Xgrid,1), size(Xgrid,2), size(offsets,1)); %will store spike counts for each presentation in a cell array
    
    thiscluster=clusters(uu,1);
    cluster_idx = find(spike_clusters==thiscluster);
    unit_spikes_all = S.kilosort(imec_num+1).spike_times(cluster_idx)/(30); %all spike times for this unit converted to milliseconds
    rf_endtime = S.rfmp1.data{size(S.rfmp1.data,1),'ns5_samps'}/(30); rf_endtime = rf_endtime(2); %time of end of last trial (in ns5 format, converted to 
    unit_spikes_all = unit_spikes_all(unit_spikes_all<=rf_endtime); %truncate data up until the end of the rf mapping session
    
    for rr=1:size(Tstimloc,1) %rows in the location matrix
        for cc=1:size(Tstimloc,2) %columns in the location matrix
            %disp(['XY = [' num2str(Xgrid(rr,cc)) ',' num2str(Ygrid(rr,cc)) ']'])
            
            locations = Tstimloc{rr,cc};
            for ll = 1:size(locations,1) %loop over locations (aka data address of each stim presentation)
                thisloc = locations(ll,:); %location is [trialnum epochnum]
                ns5_samp = S.rfmp1.data{thisloc(1),'ns5_samps'}/(30); %[start end] of trial in ns5 clock format (converted to milliseconds)
                
                stim_on = ns5_samp(1) + T{thisloc(1),'STIM_ON'}{1}(thisloc(2)); %absolute time of stimulus on in the full recording (in milliseconds)
       
                
                for oo = 1:size(offsets,1) %loop over time bins to consider
                   spike_rng = stim_on + offsets(oo,:); % sets relative time window to count spikes in
                   spike_ct = sum(unit_spikes_all>=spike_rng(1) & unit_spikes_all<spike_rng(2));
                   
                   previous_counts = rf_map_cell{rr,cc,oo};
                   new_counts = [previous_counts; spike_ct];
                   rf_map_cell{rr,cc,oo} = new_counts;
                   
                   %disppercent(oo/size(offsets,1));
                end
                
            end
        end
    end
    
    avg_rf_map = cellfun(@nanmean, rf_map_cell);
    avg_rf_map = avg_rf_map*(1000/twin); %convert avg rf map to Hz firing rates
    
    eval(['IMEC' num2str(imec_num) '.unit(uu).avg_rf_map = avg_rf_map;']);
    eval(['IMEC' num2str(imec_num) '.unit(uu).rf_map_cell = rf_map_cell;']);
    toc
    
end









%% LAZY CHECK: loop over rows of trial matrix, just pull stimulus onsets, take stored spike times (in data table) and register to each stimulus onset, and keep spikes (from 0 to 500ms) and collate across all stimuli and trials for each neuron
% to look for time-locked responses to see if data is properly registered.
imec_str = 'spiketimes_imec1'; %define the probe data to pull spikes from
n_units = size(T{1,imec_str}{1},2); %get number N of units
kept_spikes_cell=cell(n_units,1);

for uu = 1:n_units %loop over units
    disp(['mapping unit ' num2str(uu)])
    kept_spikes=[];
    
    for rr = 1:size(T,1)
       starts = T{rr,'STIM_ON'}{1};
       spikes = round(T{rr,imec_str}{1}{1,uu});
       for ss=1:numel(starts)
           thisstart = starts(ss);
           adjusted_spikes = spikes-thisstart;
           kept_spikes = [kept_spikes adjusted_spikes(adjusted_spikes>0 & adjusted_spikes<800)];
       end
    end
    kept_spikes_cell{uu,1} = kept_spikes;
end

bins=0:10:800;
unit_histos = [];
for uu = 1:n_units
    n=histc(kept_spikes_cell{uu,1},bins);
    unit_histos(uu,:) = n;
end

figure(17);
for uu=1:n_units
    clf
    plot(bins,unit_histos,'b-'); axis square; title(['unit ' num2str(uu)]); xlabel('time bin (ms)'); ylabel('total spike count');
    pause;
end;


%% More refined check: pull spikes from S.kilosort(IMEC+1) and use S.rfmp1.data{trial,'ns5_samps'} to register the stim on/off times
imec_num=1;
imec_str = ['spiketimes_imec' num2str(imec_num)]; %define the probe data to pull spikes from
n_units = size(T{1,imec_str}{1},2); %get number N of units
clusters = S.kilosort(imec_num+1).clusters{:,'cluster'};
spike_clusters = S.kilosort(imec_num+1).spike_clusters;
spike_times = double(S.kilosort(imec_num+1).spike_times/(30)); %convert spike times to milliseconds (from clock times)

for uu = 1:n_units %loop over units
    disp(['mapping unit ' num2str(uu)])
    idx = find(spike_clusters==clusters(uu,1));
    unit_spikes_all = spike_times(idx); %all spikes for this neuron
    
    kept_spikes=[];
    
    for rr = 1:size(T,1)
       ns5_times = double(T{rr,'ns5_samps'}/(30)); %start/stop times for this trial in milliseconds
        
       starts = double(T{rr,'STIM_ON'}{1}); %these are already in milliseconds
       
       for ss=1:numel(starts)
           thisstart = starts(ss) + ns5_times(1); %trial start time plus offset relative stimulus start time within this trial
           adjusted_spikes = unit_spikes_all-thisstart;
           kept_spikes = [kept_spikes adjusted_spikes(adjusted_spikes>0 & adjusted_spikes<800)];
       end
    end
    kept_spikes_cell{uu,1} = kept_spikes;
end

bins=0:10:800;
unit_histos = [];
for uu = 1:n_units
    n=histc(kept_spikes_cell{uu,1},bins);
    unit_histos(uu,:) = n;
end

figure(17);
for uu=1:n_units
    clf
    plot(bins,unit_histos(uu,:),'b-'); axis square; title(['unit ' num2str(uu)]); xlabel('time bin (ms)'); ylabel('total spike count');
    pause;
end;

figure(17);
for baseuu=[1:25:715]
    clf
    for pp=1:25
        subplot(5,5,pp);
        plot(bins,unit_histos(baseuu+(pp-1),:),'b-'); axis square; title(num2str(baseuu+(pp-1)));
    end;
    figureFormatCAH; pause;
end;




%% Simple plotting code

%for plotting one rf map
figure(13);
uu = 372;
clf
thismap = IMEC1.unit(uu).avg_rf_map;
clim = [0 max([thismap(:); 1])];

for oo = 1:25
    subplot(5,5,oo);
    imagesc(thismap(:,:,oo), clim); axis off;
    if oo==1
        title({['unit ' num2str(uu)];['peakFR ' num2str(clim(2))];[num2str(offsets(oo,:))]});
    else
        title(num2str(offsets(oo,:)));
    end
    colormap('gray');
end



%for looping over rf maps and plotting
figure(13);
for uu = 1:715
    clf
    thismap = IMEC0.unit(uu).avg_rf_map;
    clim = [0 max([thismap(:); 1])];
    
    for oo = 1:25
       subplot(5,5,oo); 
       imagesc(thismap(:,:,oo), clim); axis off; axis square;
       if oo==1
          title({['unit ' num2str(uu)];['peakFR ' num2str(clim(2))];[num2str(offsets(oo,:))]}); 
       else
          title(num2str(offsets(oo,:)));
       end
       colormap('gray');
    end
    pause(0.25);
end



%for looping over rf maps and plotting each row as a neuron (8 rows)
figure(13);
for uu = 1:8:708
    clf
    for ee = 0:7
        thismap = IMEC1.unit(uu+ee).avg_rf_map;
        clim = [0 max([thismap(:); 1])];
        for cc=1:5
            subplot(8,5,(ee*5 + cc));
            imagesc(thismap(:,:, 5*(cc-1) + 1), clim); axis off; axis square;
            if cc==1
                title({['unit ' num2str(uu+ee)];['peakFR ' num2str(clim(2))];[num2str(offsets(5*(cc-1) + 1,:))]}); 
            else
                title([num2str(offsets(5*(cc-1) + 1,:))]);
            end
        end
    end
    colormap('gray');
    pause
end

