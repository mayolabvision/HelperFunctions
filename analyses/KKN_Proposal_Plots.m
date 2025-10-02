%% --------------------------------- KKN PROPOSAL PLOTS --------------------------------- 
clear; clc;
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions'))

load('/Users/kendranoneman/Data/sapu_dualhemi/concat_tables.mat', 'C', 'Tmdir', 'Tpurs', 'Trfmp');

red_pursuity = [245,29,37];
purple_saccadey = [162,29,245];
green_visual = [24,213,115];

%% Creating concatenated tables of all sessions
data_path = '/Users/kendranoneman/Data/sapu_dualhemi';
files = dir(fullfile(data_path, '*.mat'));

S_all = cell(numel(files),4);
for i = 1:numel(files)
    fprintf('%d of %d: %s\n', i, numel(files), files(i).name);

    load(fullfile(data_path, files(i).name), 'S');
    fnames = fieldnames(S);

    % -------------------------------------- MDIR -------------------------------------- 
    tbl1 = table();
    s_fields = fieldnames(S);
    for f = 1:numel(s_fields)
        fname = s_fields{f};
        if contains(fname, {'dirmem', 'mdir'})
            this_struct = S.(fname);
            if isstruct(this_struct) && isfield(this_struct, 'tbl') && istable(this_struct.tbl)
                this_tbl = this_struct.tbl;
                if ismember('STIM8_ON', this_tbl.Properties.VariableNames)
                    this_tbl.STIM8_ON = [];
                end
                tbl1 = [tbl1; this_tbl];
            end
        end
    end

    tbl1.monkey = categorical(repmat(missing, height(tbl1), 1), {'walter','scrappy'});
    tbl1.monkey(startsWith(string(tbl1.sess_name), 'Y')) = categorical("walter");
    tbl1.monkey(startsWith(string(tbl1.sess_name), 'k')) = categorical("scrappy");

    % Only include correct trials
    tbl1 = tbl1(tbl1.result=='CORRECT',:);

    % Calculate saccade latency & toss too fast or slow trials
    tbl1.saccLatency = tbl1.SACCADE-cellfun(@(q) q(1), tbl1.FIX_OFF);
    tbl1 = tbl1(tbl1.saccLatency >= 50 & tbl1.saccLatency < 300,:);

    % Calculate eye traces & toss trials with bad first saccade
    eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x,'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',84,'PLOT_TRIAL',false), tbl1.eyedata, 'uni', 0);
    eyeVel = cellfun(@(q) calcDerivative_eyeTraces(q), cellfun(@(x) filterEyeTraces_EyeLink(x,'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',40,'PLOT_TRIAL',false), tbl1.eyedata, 'uni', 0), 'uni', 0);                                                                        
    eyeAcc = cellfun(@(q) calcDerivative_eyeTraces(q), eyeVel, 'uni', 0);
    tbl1.eyePos = eyePos; tbl1.eyeVel = eyeVel; tbl1.eyeAcc = eyeAcc;

    [th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), tbl1.eyePos, 'uni', 0);
    tbl1 = tbl1(cellfun(@(r,t,s,a) (max(r(s-10:s+75)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s-10:s+75))) - a)) < 90), rh, th, num2cell(tbl1.SACCADE), num2cell(tbl1.angle), 'uni', 1),:);
    
    frs_perTrial = (cellfun(@(u) mean(cellfun(@(q) numel(q), u, 'uni', 1)), tbl1.spiketimes_1, 'uni', 1) + cellfun(@(u) mean(cellfun(@(q) numel(q), u, 'uni', 1)), tbl1.spiketimes_2, 'uni', 1))./2;
    tbl1 = tbl1(frs_perTrial < (mean(frs_perTrial) + 3*std(frs_perTrial)) & frs_perTrial > (mean(frs_perTrial) - 3*std(frs_perTrial)),:);

    S1.mdir.tbl = tbl1;

    % -------------------------------------- PURS -------------------------------------- 
    tbl1 = table();
    s_fields = fieldnames(S);
    for f = 1:numel(s_fields)
        fname = s_fields{f};
        if contains(fname, {'purs', 'pursuit'})
            this_struct = S.(fname);
            if isstruct(this_struct) && isfield(this_struct, 'tbl') && istable(this_struct.tbl)
                tbl1 = [tbl1; this_struct.tbl];
            end
        end
    end

    tbl1.monkey = categorical(repmat(missing, height(tbl1), 1), {'walter','scrappy'});
    tbl1.monkey(startsWith(string(tbl1.sess_name), 'Y')) = categorical("walter");
    tbl1.monkey(startsWith(string(tbl1.sess_name), 'k')) = categorical("scrappy");

    % Only include correct trials
    tbl1 = tbl1(tbl1.result=='CORRECT',:);

    % Calculate eye traces 
    eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x,'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',84,'PLOT_TRIAL',false), tbl1.eyedata, 'uni', 0);
    eyeVel = cellfun(@(q) calcDerivative_eyeTraces(q), cellfun(@(x) filterEyeTraces_EyeLink(x,'SAMPLING_FREQUENCY',1000,'CUTOFF_FREQUENCY',40,'PLOT_TRIAL',false), tbl1.eyedata, 'uni', 0), 'uni', 0);                                                                        
    eyeAcc = cellfun(@(q) calcDerivative_eyeTraces(q), eyeVel, 'uni', 0);
    tbl1.eyePos = eyePos; tbl1.eyeVel = eyeVel; tbl1.eyeAcc = eyeAcc;

    % Calculate pursuit onset and pursuit type
    [pursuitOnsets,rxnTimes,msOffsets,csOnsets,csVelocities,csPeaks,csOffsets,csAngles,crossingTimes] = deal(nan(height(tbl1), 1));
    csTypes = cell(height(tbl1),1);
    for t = 1:height(tbl1)
        [pursuit_onset,rxnTime,msOffset,csOnset,csVelocity,csPeak,csOffset,csAngle,csType] = detect_pursuitOnset(tbl1.eyePos{t},tbl1.eyeVel{t},tbl1.PURSUIT_TARG_ON(t),tbl1(t,:).params.block.crossingTime,tbl1.pursuitSpeed(t),tbl1.angle(t),'PLOT_TRACES',false);
        pursuitOnsets(t) = pursuit_onset; rxnTimes(t) = rxnTime; msOffsets(t) = msOffset; csOnsets(t) = csOnset; csVelocities(t) = csVelocity; csPeaks(t) = csPeak; csOffsets(t) = csOffset; csAngles(t) = csAngle; csTypes{t} = csType;
        crossingTimes(t) = tbl1(t,:).params.block.crossingTime;
    end

    tbl1.pursuitOnset = pursuitOnsets; tbl1.pursuitLatency = rxnTimes;
    tbl1.msOffset = msOffsets; tbl1.CROSSING_TIME = crossingTimes;
    tbl1.csTimes = [csOnsets, csPeaks, csOffsets]; tbl1.csVelocity = csVelocities; tbl1.csAngle = csAngles;
    tbl1.pursType = csTypes; tbl1.pursType = categorical(string(tbl1.pursType));

    tbl1 = movevars(tbl1,{'pursuitOnset','pursuitLatency','msOffset','pursType','csTimes','csVelocity','csAngle'},'Before','result');
    tbl1 = movevars(tbl1,{'CROSSING_TIME'},'After','PURSUIT_TARG_ON');

    frs_perTrial = (cellfun(@(u) mean(cellfun(@(q) numel(q), u, 'uni', 1)), tbl1.spiketimes_1, 'uni', 1) + cellfun(@(u) mean(cellfun(@(q) numel(q), u, 'uni', 1)), tbl1.spiketimes_2, 'uni', 1))./2;
    tbl1 = tbl1(frs_perTrial < (mean(frs_perTrial) + 3*std(frs_perTrial)) & frs_perTrial > (mean(frs_perTrial) - 3*std(frs_perTrial)),:);

    S1.purs.tbl = tbl1;

    % -------------------------------------- RFMP -------------------------------------- 
    tbl1 = table();
    s_fields = fieldnames(S);
    for f = 1:numel(s_fields)
        fname = s_fields{f};
        if contains(fname, {'rfmp', 'rfMapping'})
            this_struct = S.(fname);
            if isstruct(this_struct) && isfield(this_struct, 'tbl') && istable(this_struct.tbl)
                this_tbl = this_struct.tbl;
                if ismember('STIM8_ON', this_tbl.Properties.VariableNames)
                    this_tbl.STIM8_ON = [];
                end
                tbl1 = [tbl1; this_tbl];
            end
        end
    end

    tbl1.monkey = categorical(repmat(missing, height(tbl1), 1), {'walter','scrappy'});
    tbl1.monkey(startsWith(string(tbl1.sess_name), 'Y')) = categorical("walter");
    tbl1.monkey(startsWith(string(tbl1.sess_name), 'k')) = categorical("scrappy");

    tbl1.STIM_ON(tbl1.result~='CORRECT') = cellfun(@(q) q(1:end-1), tbl1.STIM_ON(tbl1.result~='CORRECT'), 'uni', 0);
    tbl1.conditions(tbl1.result~="CORRECT") = cellfun(@(q) q(1:end-1), tbl1.conditions(tbl1.result~='CORRECT'), 'uni', 0);
    tbl1 = tbl1(~cellfun(@(q) any(isnan(q)), tbl1.STIM_OFF, 'uni', 1),:);

    S1.rfmp.tbl = tbl1;

    % ------------------------------------ CLUSTERS ------------------------------------ 
    S1.kilosort = S.kilosort;
    
    if height(S1.mdir.tbl) > 80 && height(S1.purs.tbl) > 80
        Snew = calculate_metrics_neuropixels(S1);

        commonVars = Snew.kilosort(1).clusters.Properties.VariableNames(ismember(Snew.kilosort(1).clusters.Properties.VariableNames, Snew.kilosort(2).clusters.Properties.VariableNames));
        
        K1 = Snew.kilosort(1).clusters(:, commonVars);
        K1.cluster_row = K1.cluster_id;

        K2 = Snew.kilosort(2).clusters(:, commonVars);
        K2.cluster_row = (K1.cluster_row(end) + K2.cluster_id + 1);

        Tclust = [K1; K2];

        Tclust.monkey = categorical(repmat(missing, height(Tclust), 1), {'walter','scrappy'});
        Tclust.monkey(startsWith(string(Tclust.sess_name), 'Y')) = categorical("walter");
        Tclust.monkey(startsWith(string(Tclust.sess_name), 'k')) = categorical("scrappy");

        S_all{i,1} = Tclust;
        S_all{i,2} = Snew.mdir.tbl;
        S_all{i,3} = Snew.purs.tbl;
        S_all{i,4} = Snew.rfmp.tbl;
    end
end

% Concatenate tables

% CLUSTERS
allTables = S_all(:,1);                  
allTables = allTables(~cellfun('isempty', allTables));
commonVars = allTables{1}.Properties.VariableNames;
for i = 2:numel(allTables)
    commonVars = intersect(commonVars, allTables{i}.Properties.VariableNames, 'stable');
end
for i = 1:numel(allTables)
    allTables{i} = allTables{i}(:, commonVars);
end

C = vertcat(allTables{:});

sessStr = string(C.sess_name);
C.monkey = categorical(repmat(missing, height(C), 1), {'walter','scrappy'});
C.monkey(startsWith(sessStr, 'Y')) = categorical("walter");
C.monkey(startsWith(sessStr, 'k')) = categorical("scrappy");

% MDIR
allTables = S_all(:,2);                  
allTables = allTables(~cellfun('isempty', allTables));
commonVars = allTables{1}.Properties.VariableNames;
for i = 2:numel(allTables)
    commonVars = intersect(commonVars, allTables{i}.Properties.VariableNames, 'stable');
end
for i = 1:numel(allTables)
    allTables{i} = allTables{i}(:, commonVars);
end

Tmdir = vertcat(allTables{:});

% PURS
allTables = S_all(:,3);                  
allTables = allTables(~cellfun('isempty', allTables));
commonVars = allTables{1}.Properties.VariableNames;
for i = 2:numel(allTables)
    commonVars = intersect(commonVars, allTables{i}.Properties.VariableNames, 'stable');
end
for i = 1:numel(allTables)
    allTables{i} = allTables{i}(:, commonVars);
end

Tpurs = vertcat(allTables{:});

% RFMP
allTables = S_all(:,4);                  
allTables = allTables(~cellfun('isempty', allTables));
commonVars = allTables{1}.Properties.VariableNames;
for i = 2:numel(allTables)
    commonVars = intersect(commonVars, allTables{i}.Properties.VariableNames, 'stable');
end
for i = 1:numel(allTables)
    allTables{i} = allTables{i}(:, commonVars);
end

Trfmp = vertcat(allTables{:});

Tmdir.sess_name = categorical(regexprep(string(Tmdir.sess_name),'_g0$',''));
Tpurs.sess_name = categorical(regexprep(string(Tpurs.sess_name),'_g0$',''));
Trfmp.sess_name = categorical(regexprep(string(Trfmp.sess_name),'_g0$',''));
C.sess_name = categorical(regexprep(string(C.sess_name),'_g0$',''));

save('/Users/kendranoneman/Data/sapu_dualhemi/concat_tables.mat', ...
     'C', 'Tmdir', 'Tpurs', 'Trfmp');

%% --------------------------------------------------------------------------------------

%% APPLYING THRESHOLDS TO CLUSTERS BEFORE MOVING FORWARD WITH PLOTTING

% Removing any units with NaN or Inf in selectivity metrics
C1 = C(~isnan(C.VMI) & ~isinf(C.VMI) & ~isnan(C.VMIdp) & ~isinf(C.VMIdp) & ~isnan(C.SPI) & ...
       ~isinf(C.SPI) & ~isnan(C.SPIdp) & ~isinf(C.SPIdp) & ~isnan(C.SPI) & ~isinf(C.SPI) & ...
       C.vis_sel_dir>0 & C.vis_sel_dir<1 & C.sac_sel_dir>0 & C.sac_sel_dir<1 & C.pur_sel_dir>0 & C.pur_sel_dir<1,:);

% Removing any units that don't fire at least once on most of the trials
PR = 0.95;
C2 = C1(C1.ratio_mdir_trials >= PR & C1.ratio_purs_trials >=PR,:);

% Minimum dp between evoked and spontaneous activity
DP = 0.3;
C3 = C2(C2.dp_vis >= DP | C2.dp_sac >= DP | C2.dp_pur >= DP,:);

% Minimum FR
FR = 1;
C4 = C3(C3.mdir_delayFR_peakDirFR >= FR | C3.purs_targFR_peakDirFR >= FR,:);

% Minimum SNR
SNR = 2;
C5 = C4(C4.snr>=SNR,:);

CC = C5;
clustCounts = groupcounts(CC, {'monkey'});

%% CHAPTER 3 - SACCADE-PURSUIT INTERACTIONS IN ONE HEMI

%% PSTH colored by VMI and SPI quantiles
nBins = 5;
sigma = 10;
visWindow = [-100,300];
sacWindow = [-200,200];
purWindow = [-200,200];

CC.VMI_quantile = discretize(CC.VMI,  quantile(CC.VMI, linspace(0, 1, nBins+1)));
CC.VMI_quantile(CC.VMI == max(CC.VMI)) = nBins;

[vis_mn_vmi,vis_sem_vmi,sac_mn_vmi,sac_sem_vmi,pur_mn_vmi,pur_sem_vmi] = deal(cell(numel(these_sess),nBins));
for q = 1:nBins
    fprintf('Bin %d of %d...\n', q, nBins);
    these_units = CC(CC.VMI_quantile==q,:);

    these_sess = unique(these_units.sess_name);
    for s = 1:numel(these_sess)
        this_mdirTbl = Tmdir(Tmdir.sess_name==these_sess(s),:);
        this_pursTbl = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
        
        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_mdirTbl.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_mdirTbl.spiketimes_2, 'uni', 0);

        % VISUAL PSTH
        time_window = visWindow;
        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes1, this_mdirTbl.TARG_ON, 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes2, this_mdirTbl.TARG_ON, 'uni', 0);

        spike_times = [vertcat(spike_times1{:}) vertcat(spike_times2{:})];
        tstep = 1;
        time = tstep + time_window(1) : tstep : time_window(2);
        
        [vis_mn, vis_sem] = deal(zeros(size(spike_times,2), length(time)));
        for unit = 1:size(spike_times,2)
            sdf = zeros(size(spike_times,1), length(time));
            for iTrial = 1:size(spike_times,1)
                spks = spike_times{iTrial,unit};
                if size(spks,2)==1
                    spks = spks';
                end
                if isempty(spks)
                    continue;
                end
                
                gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
                sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
            end

           [mn, sem, ~, ~] = sem_errorbar(sdf .* 1000);
           vis_mn(unit,:) = mn;
           vis_sem(unit,:) = sem;
        end

        vis_mn_vmi{s,q} = vis_mn;
        vis_sem_vmi{s,q} = vis_sem;

        % SAC PSTH
        time_window = sacWindow;
        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_mdirTbl.SACCADE), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_mdirTbl.SACCADE), 'uni', 0);

        spike_times = [vertcat(spike_times1{:}) vertcat(spike_times2{:})];
        tstep = 1;
        time = tstep + time_window(1) : tstep : time_window(2);
        
        [sac_mn, sac_sem] = deal(zeros(size(spike_times,2), length(time)));
        for unit = 1:size(spike_times,2)
            sdf = zeros(size(spike_times,1), length(time));
            for iTrial = 1:size(spike_times,1)
                spks = spike_times{iTrial,unit};
                if size(spks,2)==1
                    spks = spks';
                end
                if isempty(spks)
                    continue;
                end
                
                gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
                sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
            end

           [mn, sem, ~, ~] = sem_errorbar(sdf .* 1000);
           sac_mn(unit,:) = mn;
           sac_sem(unit,:) = sem;
        end

        sac_mn_vmi{s,q} = sac_mn;
        sac_sem_vmi{s,q} = sac_sem;

        % PURS PSTH
        time_window = purWindow;

        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl.spiketimes_2, 'uni', 0);

        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);

        spike_times = [vertcat(spike_times1{:}) vertcat(spike_times2{:})];
        tstep = 1;
        time = tstep + time_window(1) : tstep : time_window(2);
        
        [pur_mn, pur_sem] = deal(zeros(size(spike_times,2), length(time)));
        for unit = 1:size(spike_times,2)
            sdf = zeros(size(spike_times,1), length(time));
            for iTrial = 1:size(spike_times,1)
                spks = spike_times{iTrial,unit};
                if size(spks,2)==1
                    spks = spks';
                end
                if isempty(spks)
                    continue;
                end
                
                gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
                sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
            end

           [mn, sem, ~, ~] = sem_errorbar(sdf .* 1000);
           pur_mn(unit,:) = mn;
           pur_sem(unit,:) = sem;
        end

        pur_mn_vmi{s,q} = pur_mn;
        pur_sem_vmi{s,q} = pur_sem;

    end
end

CC.SPI_quantile = discretize(CC.SPI, quantile(CC.SPI, linspace(0, 1, nBins+1)));
CC.SPI_quantile(CC.SPI == max(CC.SPI)) = nBins;

[vis_mn_spi,vis_sem_spi,sac_mn_spi,sac_sem_spi,pur_mn_spi,pur_sem_spi] = deal(cell(numel(these_sess),nBins));
for q = 1:nBins
    fprintf('Bin %d of %d...\n', q, nBins);
    these_units = CC(CC.SPI_quantile==q,:);

    these_sess = unique(these_units.sess_name);
    for s = 1:numel(these_sess)
        this_mdirTbl = Tmdir(Tmdir.sess_name==these_sess(s),:);
        this_pursTbl = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
        
        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_mdirTbl.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_mdirTbl.spiketimes_2, 'uni', 0);

        % VISUAL PSTH
        time_window = visWindow;
        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes1, this_mdirTbl.TARG_ON, 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes2, this_mdirTbl.TARG_ON, 'uni', 0);

        spike_times = [vertcat(spike_times1{:}) vertcat(spike_times2{:})];
        tstep = 1;
        time = tstep + time_window(1) : tstep : time_window(2);
        
        [vis_mn, vis_sem] = deal(zeros(size(spike_times,2), length(time)));
        for unit = 1:size(spike_times,2)
            sdf = zeros(size(spike_times,1), length(time));
            for iTrial = 1:size(spike_times,1)
                spks = spike_times{iTrial,unit};
                if size(spks,2)==1
                    spks = spks';
                end
                if isempty(spks)
                    continue;
                end
                
                gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
                sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
            end

           [mn, sem, ~, ~] = sem_errorbar(sdf .* 1000);
           vis_mn(unit,:) = mn;
           vis_sem(unit,:) = sem;
        end

        vis_mn_spi{s,q} = vis_mn;
        vis_sem_spi{s,q} = vis_sem;

        % SAC PSTH
        time_window = sacWindow;
        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_mdirTbl.SACCADE), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_mdirTbl.SACCADE), 'uni', 0);

        spike_times = [vertcat(spike_times1{:}) vertcat(spike_times2{:})];
        tstep = 1;
        time = tstep + time_window(1) : tstep : time_window(2);
        
        [sac_mn, sac_sem] = deal(zeros(size(spike_times,2), length(time)));
        for unit = 1:size(spike_times,2)
            sdf = zeros(size(spike_times,1), length(time));
            for iTrial = 1:size(spike_times,1)
                spks = spike_times{iTrial,unit};
                if size(spks,2)==1
                    spks = spks';
                end
                if isempty(spks)
                    continue;
                end
                
                gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
                sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
            end

           [mn, sem, ~, ~] = sem_errorbar(sdf .* 1000);
           sac_mn(unit,:) = mn;
           sac_sem(unit,:) = sem;
        end

        sac_mn_spi{s,q} = sac_mn;
        sac_sem_spi{s,q} = sac_sem;

        % PURS PSTH
        time_window = purWindow;

        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl.spiketimes_2, 'uni', 0);

        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);

        spike_times = [vertcat(spike_times1{:}) vertcat(spike_times2{:})];
        tstep = 1;
        time = tstep + time_window(1) : tstep : time_window(2);
        
        [pur_mn, pur_sem] = deal(zeros(size(spike_times,2), length(time)));
        for unit = 1:size(spike_times,2)
            sdf = zeros(size(spike_times,1), length(time));
            for iTrial = 1:size(spike_times,1)
                spks = spike_times{iTrial,unit};
                if size(spks,2)==1
                    spks = spks';
                end
                if isempty(spks)
                    continue;
                end
                
                gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
                sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
            end

           [mn, sem, ~, ~] = sem_errorbar(sdf .* 1000);
           pur_mn(unit,:) = mn;
           pur_sem(unit,:) = sem;
        end

        pur_mn_spi{s,q} = pur_mn;
        pur_sem_spi{s,q} = pur_sem;
    end
end


%%
spi_cmap = customColormap(red_pursuity,purple_saccadey, nBins);
vmi_cmap = customColormap(purple_saccadey, green_visual, nBins);

[psth_mn,psth_sem] = deal(zeros(nBins,1200));
for q = 1:nBins
    psth_mn_all = [vertcat(vis_mn_vmi{:,q}) vertcat(sac_mn_vmi{:,q}) vertcat(pur_mn_vmi{:,q})];
    % psth_mn_norm = psth_mn_all ./ max(psth_mn_all, [], 2);
    % 
    % row_min = min(psth_mn_all, [], 2);
    % row_max = max(psth_mn_all, [], 2);
    % psth_mn_norm = (psth_mn_all - row_min) ./ (row_max - row_min);

    psth_mn_norm = zscore(psth_mn_all,[],2);

    psth_mn(q,:) = mean(psth_mn_norm,1);
    
    psth_sem_all = [vertcat(vis_sem_vmi{:,q}) vertcat(sac_sem_vmi{:,q}) vertcat(pur_sem_vmi{:,q})];
    
    % psth_sem_norm = psth_sem_all ./ max(psth_mn_all, [], 2);
    % psth_sem_norm =  psth_sem_all ./ (row_max - row_min);

    psth_sem_norm = psth_sem_all ./ std(psth_mn_all, 0, 2);  % divide by the same std as z-scoring

    psth_sem(q,:) = mean(psth_sem_norm,1);
end

f3a = figure;
f3a.Position = [100 100 1000 400];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

ax1(1) = nexttile;
x = visWindow(1):(visWindow(2)-1);

%xline(0,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold on;
%fill([[50,150] fliplr([50,150])], [[1,1] fliplr([0,0])], [0.7,0.7,0.7], 'linestyle', 'none', 'FaceAlpha', 0.1);
for q = 1:nBins
    y = psth_mn(q,1:400);
    yu = psth_mn(q,1:400)+psth_sem(q,1:400); yl = psth_mn(q,1:400)-psth_sem(q,1:400);

    fill([x fliplr(x)], [yu fliplr(yl)], vmi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',vmi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

ax1(2) = nexttile;
x = sacWindow(1):(sacWindow(2)-1);

%xline(0,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold on;
%fill([[-50,50] fliplr([-50,50])], [[1,1] fliplr([0,0])], [0.7,0.7,0.7], 'linestyle', 'none', 'FaceAlpha', 0.1);
for q = 1:nBins
    y = psth_mn(q,401:800);
    yu = psth_mn(q,401:800)+psth_sem(q,401:800); yl = psth_mn(q,401:800)-psth_sem(q,401:800);

    fill([x fliplr(x)], [yu fliplr(yl)], vmi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',vmi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

nexttile(3);
hold on

edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
% Create dummy lines for the legend (they won't be plotted in the main figure)
for q = 1:nBins
    plot(nan, nan, '-', 'Color', vmi_cmap(q,:), 'LineWidth', 3, 'DisplayName', sprintf('VMI = [%.2f, %.2f]', edges(q), edges(q+1)));
end

lgd = legend('Location', 'best', 'Orientation', 'vertical');
lgd.FontSize = 14;
lgd.ItemTokenSize = [30, 18];
axis off  % hide axes since this tile is only for the legend


ax1(3) = nexttile(4);
x = purWindow(1):(purWindow(2)-1);

%xline(0,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold on;
%fill([[-50,50] fliplr([-50,50])], [[1,1] fliplr([0,0])], [0.7,0.7,0.7], 'linestyle', 'none', 'FaceAlpha', 0.1);
for q = 1:nBins
    y = psth_mn(q,801:end);
    yu = psth_mn(q,801:end)+psth_sem(q,801:end); yl = psth_mn(q,801:end)-psth_sem(q,801:end);

    fill([x fliplr(x)], [yu fliplr(yl)], vmi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',vmi_cmap(q,:),'LineWidth',3);
end
prettyFig;

linkaxes(ax1,'y')

savebigPDF(f3a, '/Users/kendranoneman/Milestones/proposal/figs/psth_vmi.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[psth_mn,psth_sem] = deal(zeros(nBins,1200));
for q = 1:nBins
    psth_mn_all = [vertcat(vis_mn_spi{:,q}) vertcat(sac_mn_spi{:,q}) vertcat(pur_mn_spi{:,q})];
    % psth_mn_norm = psth_mn_all ./ max(psth_mn_all, [], 2);
    % 
    %row_min = min(psth_mn_all, [], 2);
    %row_max = max(psth_mn_all, [], 2);
    %psth_mn_norm = (psth_mn_all - row_min) ./ (row_max - row_min);

    psth_mn_norm = zscore(psth_mn_all,[],2);

    psth_mn(q,:) = mean(psth_mn_norm,1);
    
    psth_sem_all = [vertcat(vis_sem_spi{:,q}) vertcat(sac_sem_spi{:,q}) vertcat(pur_sem_spi{:,q})];
    
    % psth_sem_norm = psth_sem_all ./ max(psth_mn_all, [], 2);
    %psth_sem_norm =  psth_sem_all ./ (row_max - row_min);

    psth_sem_norm = psth_sem_all ./ std(psth_mn_all, 0, 2);  % divide by the same std as z-scoring

    psth_sem(q,:) = mean(psth_sem_norm,1);
end

f3b = figure;
f3b.Position = [100 100 1000 400];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

ax1(1) = nexttile;
x = visWindow(1):(visWindow(2)-1);

%xline(0,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold on;
%fill([[50,150] fliplr([50,150])], [[1,1] fliplr([0,0])], [0.7,0.7,0.7], 'linestyle', 'none', 'FaceAlpha', 0.1);
for q = 1:nBins
    y = psth_mn(q,1:400);
    yu = psth_mn(q,1:400)+psth_sem(q,1:400); yl = psth_mn(q,1:400)-psth_sem(q,1:400);

    fill([x fliplr(x)], [yu fliplr(yl)], spi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',spi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

ax1(2) = nexttile;
x = sacWindow(1):(sacWindow(2)-1);

%xline(0,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold on;
%fill([[-50,50] fliplr([-50,50])], [[1,1] fliplr([0,0])], [0.7,0.7,0.7], 'linestyle', 'none', 'FaceAlpha', 0.1);
for q = 1:nBins
    y = psth_mn(q,401:800);
    yu = psth_mn(q,401:800)+psth_sem(q,401:800); yl = psth_mn(q,401:800)-psth_sem(q,401:800);

    fill([x fliplr(x)], [yu fliplr(yl)], spi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',spi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

nexttile(3);
hold on

edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
% Create dummy lines for the legend (they won't be plotted in the main figure)
for q = 1:nBins
    plot(nan, nan, '-', 'Color', spi_cmap(q,:), 'LineWidth', 3, 'DisplayName', sprintf('SPI = [%.2f, %.2f]', edges(q), edges(q+1)));
end

lgd = legend('Location', 'best', 'Orientation', 'vertical');
lgd.FontSize = 14;
lgd.ItemTokenSize = [30, 18];
axis off  % hide axes since this tile is only for the legend


ax1(3) = nexttile(4);
x = purWindow(1):(purWindow(2)-1);

%xline(0,'--','color',[0.5,0.5,0.5],'linewidth',2)
hold on;
%fill([[-50,50] fliplr([-50,50])], [[1,1] fliplr([0,0])], [0.7,0.7,0.7], 'linestyle', 'none', 'FaceAlpha', 0.1);
for q = 1:nBins
    y = psth_mn(q,801:end);
    yu = psth_mn(q,801:end)+psth_sem(q,801:end); yl = psth_mn(q,801:end)-psth_sem(q,801:end);

    fill([x fliplr(x)], [yu fliplr(yl)], spi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',spi_cmap(q,:),'LineWidth',3);
end
prettyFig;

linkaxes(ax1,'y')
savebigPDF(f3b, '/Users/kendranoneman/Milestones/proposal/figs/psth_spi.pdf')

%% VMI and SPI as a function of depth

spi_cmap = customColormap(red_pursuity, purple_saccadey, 256);
vmi_cmap = customColormap(purple_saccadey, green_visual, 256);

f3c = figure;
f3c.Position = [100 100 500 500];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

for p = 1:2
    ax1(p) = nexttile;
    
    if p==1
        xx = CC.VMI;
        cmap = vmi_cmap;
    else
        xx = CC.SPI;
        cmap = spi_cmap;
    end

    zz = cellfun(@(q) q(2), CC.unit_locations, 'uni', 1) ./1000;
    zz = zz - (CC.probe_depth_mm);
    
    % Scatter plot with color mapped to xx
    scatter(ax1(p), xx, zz, 10, xx, 'filled');  
    hold on;
    
    % --- Linear fit line ---
    lm = fitlm(xx, zz);  % fit linear model
    xfit = linspace(min(xx), max(xx), 100);
    yfit = predict(lm, xfit');
    plot(ax1(p), xfit, yfit, 'k-', 'LineWidth', 2);
    
    % --- Correlation stats ---
    [r, pval] = corr(xx, zz, 'Type', 'Spearman');
    
    % Add text to plot (top-left corner)
    text(min(xx), max(zz), sprintf('r = %.2f\np = %.3g', r, pval), ...
        'VerticalAlignment', 'top', 'FontSize', 10, 'Parent', ax1(p));
    
    % --- Colormap and colorbar ---
    colormap(ax1(p), cmap);
    cb = colorbar(ax1(p));
    %clim([-1,1]);
    
    prettyFig;
end

linkaxes(ax1,'xy')

%savebigPDF(f3c, '/Users/kendranoneman/Milestones/proposal/figs/vmi_spi_depth.pdf')


%% Ternary plot of FR 

mnFR = [CC.vis_meanFR CC.sac_meanFR CC.pur_meanFR];
normFR = mnFR ./ sum(mnFR,2);

% normFR is 1000x3, rows sum to 1

a = normFR(:,1); % Visual
b = normFR(:,2); % Motor
c = normFR(:,3); % Pursuit

% Barycentric -> Cartesian mapping (since a+b+c=1)
x = b + 0.5*c;
y = (sqrt(3)/2)*c;

% Define number of bins
nbins = 50;

% Fix the bin edges to cover the full ternary triangle
xedges = linspace(0, 1, nbins+1);
yedges = linspace(0, sqrt(3)/2, nbins+1);

% 2D histogram
[counts, Xedges, Yedges] = histcounts2(x, y, xedges, yedges);

% Run one-sample t-tests against 0.5
[~,pA,~,statsA] = ttest(a, 0.3333333333);
[~,pB,~,statsB] = ttest(b, 0.3333333333);
[~,pC,~,statsC] = ttest(c, 0.3333333333);

fprintf('Visual: mean = %.3f, t(%d)=%.2f, p=%.3g\n', mean(a), statsA.df, statsA.tstat, pA);
fprintf('Motor:  mean = %.3f, t(%d)=%.2f, p=%.3g\n', mean(b), statsB.df, statsB.tstat, pB);
fprintf('Pursuit: mean = %.3f, t(%d)=%.2f, p=%.3g\n', mean(c), statsC.df, statsC.tstat, pC);

% Plot
f3d = figure;
imagesc(Xedges, Yedges, counts'); 
axis xy equal tight
colormap gray
colorbar
title('Ternary heatmap of unit distribution');

% Overlay triangle edges
hold on;
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2);

% Corner labels
text(0,0,'Visual','HorizontalAlignment','center','VerticalAlignment','top','FontSize',12)
text(1,0,'Motor','HorizontalAlignment','center','VerticalAlignment','top','FontSize',12)
text(0.5,sqrt(3)/2,'Pursuit','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12)

prettyFig;
savebigPDF(f3d, '/Users/kendranoneman/Milestones/proposal/figs/fr_ternary_heatmap.pdf')

%% Distributions of SPI and VMI (per monkey)

f3e = figure;
f3e.Position = [100 100 1300 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

ax1(1) = nexttile;
values = CC.VMI(CC.monkey=='scrappy');
histStyle_KKN(values, 'BIN_WIDTH', 0.1, 'Y_LIMITS', [0 600], 'X_LIMITS', [-1,1]);

ax1(2) = nexttile;
values = CC.SPI(CC.monkey=='scrappy');
histStyle_KKN(values, 'BIN_WIDTH', 0.1, 'Y_LIMITS', [0 600], 'X_LIMITS', [-1,1]);

ax1(3) = nexttile;
values = CC.VMI(CC.monkey=='walter');
histStyle_KKN(values, 'BIN_WIDTH', 0.1, 'Y_LIMITS', [0 110], 'X_LIMITS', [-1,1]);

ax1(4) = nexttile;
values = CC.SPI(CC.monkey=='walter');
histStyle_KKN(values, 'BIN_WIDTH', 0.1, 'Y_LIMITS', [0 110], 'X_LIMITS', [-1,1]);

prettyFig;
%savebigPDF(f3e, '/Users/kendranoneman/Milestones/proposal/figs/hists_spi_vmi.pdf')

%% VMI vs. SPI

f3f = figure;
f3f.Position = [100 100 500 500];

% Extract variables
x = CC.VMI;
y = CC.SPI;

% Define edges (here I fix them to [-1,1], but you can use min/max instead)
nbins = 44;
xedges = linspace(min(x), max(x), nbins+1);
yedges = linspace(min(y), max(y), nbins+1);

% Compute 2D histogram
[counts,~,~] = histcounts2(x, y, xedges, yedges);

% Plot heatmap
imagesc(xedges, yedges, counts');  
axis xy; axis equal tight
colormap(gray);
colorbar;

hold on;

% --- Linear fit line ---
lm = fitlm(x, y);  % fit linear model
xfit = linspace(min(x), max(y), 100);
yfit = predict(lm, xfit');
plot(xfit, yfit, 'r-', 'LineWidth', 1);

% --- Correlation stats ---
[r, pval] = corr(x, y, 'Type', 'Spearman');

xlim([-1,1]);
ylim([-1,1]);

prettyFig;

savebigPDF(f3f, '/Users/kendranoneman/Milestones/proposal/figs/spi_versus_vmi.pdf')


%% RSC

sess = unique(CC.sess_name);
pairs_all = cell(length(sess),1);

for s = 1:length(sess)
    fprintf('%d of %d\n', s, numel(sess));
    this_mdir = Tmdir(Tmdir.sess_name==sess(s),:);
    this_purs = Tpurs(Tpurs.sess_name==sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);

    pairs_probes = cell(4,1);
    for p1 = 1:2
        tStart = tic;  % start timer
        for p2 = 1:2
            fprintf('%d - %d\n', p1, p2);

            p1_units = CC(CC.sess_name==sess(s) & CC.probe_index==p1,:);
            p2_units = CC(CC.sess_name==sess(s) & CC.probe_index==p2,:);
        
            pairs = cell(height(p1_units)*height(p2_units),24);
            rr = 1;
            for n1 = 1:height(p1_units)
                for n2 = 1:height(p2_units)
                    % MDIR
                    angs = sort(unique(this_mdir.angle));
                    amps = sort(unique(this_mdir.distance));

                    if 0.5*(p1_units.vis_sel_dir(n1)+p1_units.sac_sel_dir(n1)) > 0.5*(p2_units.vis_sel_dir(n2)+p2_units.sac_sel_dir(n2))
                        bestDir = p1_units.mdir_delayFR_peakDir(n1);
                    else
                        bestDir = p2_units.mdir_delayFR_peakDir(n2);
                    end

                    relAngs = mod(angs - bestDir + 180, 360) - 180;
                    [~, sortIdx] = sort(relAngs);

                    [rsc_perDir_mdir, fr_perDir_mdir] = deal(zeros(1, numel(angs)));
                    n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                    for d = 1:numel(angs)
                        thisAng = angs(d);

                        n12_fr = []; n1_zfr = []; n2_zfr = [];
                        for a = 1:numel(amps)
                            n12_fr = [n1_fr; [cellfun(@(q,w,v) ((sum(q{n1}>w(1) & q{n1}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p1))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1);cellfun(@(q,w,v) ((sum(q{n2}>w(1) & q{n2}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p2))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1)]];
                            
                            n1_zfr = [n1_zfr; zscore(cellfun(@(q,w,v) ((sum(q{n1}>w(1) & q{n1}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p1))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1))];
                            n2_zfr = [n2_zfr; zscore(cellfun(@(q,w,v) ((sum(q{n2}>w(1) & q{n2}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p2))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1))];
                        end
                        n1_zfr_all = [n1_zfr_all; n1_zfr];
                        n2_zfr_all = [n2_zfr_all; n2_zfr];
                        n12_fr_all = [n12_fr_all; n12_fr];

                        [rho,~] = corr(n1_zfr, n2_zfr);

                        rsc_perDir_mdir(sortIdx(d)) = rho;
                        fr_perDir_mdir(sortIdx(d)) = mean(n12_fr);
                    end

                    n1_fr = cellfun(@(q,w,v) ((sum(q{n1}>w(1) & q{n1}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p1)), this_mdir.TARG_ON, num2cell(this_mdir.SACCADE), 'uni', 1);
                    n2_fr = cellfun(@(q,w,v) ((sum(q{n2}>w(1) & q{n2}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p2)), this_mdir.TARG_ON, num2cell(this_mdir.SACCADE), 'uni', 1);

                    [rsig_mdir,~] = corr(n1_fr,n2_fr);

                    [rsc_mdir,~] = corr(n1_zfr_all, n1_zfr_all);
                    mnFR_mdir = mean(n12_fr_all);

                    % PURS
                    angs = sort(unique(this_purs.angle));
                    amps = sort(unique(this_purs.pursuitSpeed));

                    if p1_units.pur_sel_dir(n1) > p2_units.pur_sel_dir(n2)
                        bestDir = p1_units.purs_targFR_peakDir(n1);
                    else
                        bestDir = p2_units.purs_targFR_peakDir(n2);
                    end

                    relAngs = mod(angs - bestDir + 180, 360) - 180;
                    [relAngsSorted, sortIdx] = sort(relAngs);

                    [rsc_perDir_purs, fr_perDir_purs] = deal(zeros(1, numel(angs)));
                    n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                    for d = 1:numel(angs)
                        thisAng = angs(d);

                        n12_fr = []; n1_zfr = []; n2_zfr = [];
                        for a = 1:numel(amps)
                            n12_fr = [n1_fr; [cellfun(@(q,w) ((sum(q{n1}>w & q{n1}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p1))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1); cellfun(@(q,w) ((sum(q{n2}>w & q{n2}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p2))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1)]];
                            
                            n1_zfr = [n1_zfr; zscore(cellfun(@(q,w) ((sum(q{n1}>w & q{n1}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p1))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1))];
                            n2_zfr = [n2_zfr; zscore(cellfun(@(q,w) ((sum(q{n2}>w & q{n2}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p2))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1))];
                        end
                        n1_zfr_all = [n1_zfr_all; n1_zfr];
                        n2_zfr_all = [n2_zfr_all; n2_zfr];
                        n12_fr_all = [n12_fr_all; n12_fr];

                        [rho,~] = corr(n1_zfr, n2_zfr);

                        rsc_perDir_purs(sortIdx(d)) = rho;
                        fr_perDir_purs(sortIdx(d)) = mean(n12_fr);
                    end

                    n1_fr = cellfun(@(q,w) ((sum(q{n1}>w & q{n1}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p1)), num2cell(this_purs.PURSUIT_TARG_ON), 'uni', 1);
                    n2_fr = cellfun(@(q,w) ((sum(q{n2}>w & q{n2}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p2)), num2cell(this_purs.PURSUIT_TARG_ON), 'uni', 1);

                    [rsig_purs,~] = corr(n1_fr,n2_fr);

                    [rsc_purs,~] = corr(n1_zfr_all, n1_zfr_all);
                    mnFR_purs = mean(n12_fr_all);
                    
                    pairs(rr,:) = {p1_units.monkey(n1), p1_units.sess_name(n1) ...
                            p1_units.probe_index(n1), p2_units.probe_index(n2), ...
                            p1_units.cluster_id(n1), p2_units.cluster_id(n2), ...
                            p1_units.VMI_quantile(n1), p2_units.VMI_quantile(n2), ...
                            p1_units.SPI_quantile(n1), p2_units.SPI_quantile(n2), ...
                            abs(p1_units.unit_locations{n1}(2)-p2_units.unit_locations{n2}(2)), ...
                            abs(mod(p1_units.vis_pref_dir(n1) - p2_units.vis_pref_dir(n2) + 180, 360) - 180), ...
                            abs(mod(p1_units.sac_pref_dir(n1) - p2_units.sac_pref_dir(n2) + 180, 360) - 180), ...
                            abs(mod(p1_units.pur_pref_dir(n1) - p2_units.pur_pref_dir(n2) + 180, 360) - 180), ...
                            mnFR_mdir, rsc_mdir, rsig_mdir, rsc_perDir_mdir, fr_perDir_mdir, ...
                            mnFR_purs, rsc_purs, rsig_purs, rsc_perDir_purs, fr_perDir_purs
                            };

                    rr = rr + 1;
                end
            end   
            pairs_probes{(p1-1)*2 + p2} = pairs;
        end
    end
    pairs_all{s} = pairs_probes;
end

f1a = figure;
f1a.Position = [100 100 1500 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

ax1(1) = nexttile;
values = Cnew.SPI(Cnew.monkey=='scrappy' & Cnew.probe_index==1);
histStyle_KKN(values, 'BIN_WIDTH', 0.1);

ax1(2) = nexttile;
values = Cnew.SPI(Cnew.monkey=='scrappy' & Cnew.probe_index==2);
histStyle_KKN(values, 'BIN_WIDTH', 0.1);%;, 'Y_LIMITS', [0 3600], 'X_LIMITS', [0,10]);

ax1(3) = nexttile;
values = Cnew.SPI(Cnew.monkey=='walter' & Cnew.probe_index==1);
histStyle_KKN(values, 'BIN_WIDTH', 0.1);

ax1(4) = nexttile;
values = Cnew.SPI(Cnew.monkey=='walter' & Cnew.probe_index==2);
histStyle_KKN(values, 'BIN_WIDTH', 0.1);

%% CHAPTER 4 - SACCADE-PURSUIT INTERACTIONS ACROSS HEMIS
%% Distributions of preferred directions
bw = 15;
face_alpha = 0.1;

f1a = figure;
f1a.Position = [100 100 1300 900]; 
tl = tiledlayout(2,2); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'tight'; 

% Scrappy 
RLIM = 300;

ax1(1) = nexttile;
values = CC.vis_pref_dir(CC.monkey=='scrappy' & CC.probe_index==1);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [10,128,255]./255, 'RLIM', RLIM, ...
                   'LINE_COLOR', [10,128,255]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.sac_pref_dir(CC.monkey=='scrappy' & CC.probe_index==1);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [128,255,10]./255, ...
                  'LINE_COLOR', [128,255,10]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.pur_pref_dir(CC.monkey=='scrappy' & CC.probe_index==1);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [108,20,245]./255, ...
                  'LINE_COLOR', [108,20,245]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
prettyFig; 

ax1(2) = nexttile;
values = CC.vis_pref_dir(CC.monkey=='scrappy' & CC.probe_index==2);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [10,128,255]./255, 'RLIM', RLIM, ...
                  'LINE_COLOR', [10,128,255]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.sac_pref_dir(CC.monkey=='scrappy' & CC.probe_index==2);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [128,255,10]./255, ...
                  'LINE_COLOR', [128,255,10]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.pur_pref_dir(CC.monkey=='scrappy' & CC.probe_index==2);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [108,20,245]./255, ...
                  'LINE_COLOR', [108,20,245]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
prettyFig; 

% Walter
RLIM = 55;

ax1(3) = nexttile;
values = CC.vis_pref_dir(CC.monkey=='walter' & CC.probe_index==1);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [10,128,255]./255, 'RLIM', RLIM, ...
                  'LINE_COLOR', [10,128,255]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.sac_pref_dir(CC.monkey=='walter' & CC.probe_index==1);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [128,255,10]./255, ...
                  'LINE_COLOR', [128,255,10]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.pur_pref_dir(CC.monkey=='walter' & CC.probe_index==1);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [108,20,245]./255, ...
                  'LINE_COLOR', [108,20,245]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
prettyFig; 

ax1(4) = nexttile;
values = CC.vis_pref_dir(CC.monkey=='walter' & CC.probe_index==2);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [10,128,255]./255, 'RLIM', RLIM, ...
                  'LINE_COLOR', [10,128,255]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.sac_pref_dir(CC.monkey=='walter' & CC.probe_index==2);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [128,255,10]./255, ...
                  'LINE_COLOR', [128,255,10]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = CC.pur_pref_dir(CC.monkey=='walter' & CC.probe_index==2);
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [108,20,245]./255, ...
                  'LINE_COLOR', [108,20,245]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
prettyFig; 


%% 

bw = 15;
face_alpha = 0.1;

f1a = figure;
f1a.Position = [100 100 1300 900];  

% Scrappy 
RLIM = 600;

values = [CC.vis_pref_dir(CC.probe_index==1); wrapTo360(CC.vis_pref_dir(CC.probe_index==2)+180)];
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [10,128,255]./255, 'RLIM', RLIM, ...
                   'LINE_COLOR', [10,128,255]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = [CC.sac_pref_dir(CC.probe_index==1); wrapTo360(CC.sac_pref_dir(CC.probe_index==2)+180)];
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [128,255,10]./255, ...
                  'LINE_COLOR', [128,255,10]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = [CC.pur_pref_dir(CC.probe_index==1); wrapTo360(CC.pur_pref_dir(CC.probe_index==2)+180)];
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', [108,20,245]./255, ...
                  'LINE_COLOR', [108,20,245]./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
prettyFig; 



%% OTHER PLOTS NOT USING RN

%% 3D - VMI and SPI as a function of depth??
spi_cmap = customColormap(red_pursuity,purple_saccadey, 256);
vmi_cmap = customColormap(purple_saccadey, green_visual, 256);

f3c = figure;
f3c.Position = [100 100 1000 400];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

for i = 1:2
    ax1(i) = nexttile;
    monkey = 'scrappy';
    probe = i;
    
    cc = CC.VMI(CC.monkey==monkey & CC.probe_index==probe);
    zz = cellfun(@(q) q(2), CC.unit_locations(CC.monkey==monkey & CC.probe_index==probe), 'uni', 1);
    zz = zz - (CC.probe_depth_mm(CC.monkey==monkey & CC.probe_index==probe) * 1000);
    
    xy = CC.probe_gridHole(CC.monkey==monkey & CC.probe_index==probe); 
    
    N = numel(xy);
    out = zeros(N,2);
    for i = 1:N
        vals = xy{i};  % 2x1 cell array of strings
    
        % --- reorder so ML/m/l is first, AP/p/a is second
        if contains(vals{1}, {'AP','p','a'}, 'IgnoreCase', true)
            vals = flipud(vals); % swap order
        end
    
        % --- convert each string to number
        for j = 1:2
            s = vals{j};
    
            if strcmpi(s,'ML') || strcmpi(s,'AP')
                num = 0;
            else
                prefix = lower(s(1));
                val = str2double(s(2:end));
                if prefix == 'm'
                    num = val;
                elseif prefix == 'l'
                    num = -val;
                elseif prefix == 'a'
                    num = val;
                elseif prefix == 'p'
                    num = -val;
                else
                    num = NaN; % fallback if unexpected string
                end
            end
            out(i,j) = num;
        end
    end
    
    xx = out(:,2);
    yy = out(:,1);
    
    % Define grid to interpolate onto
    xlin = linspace(min(xx), max(xx), 50);   % 50 points along X
    ylin = linspace(min(yy), max(yy), 50);   % 50 points along Y
    [Xgrid, Ygrid] = meshgrid(xlin, ylin);
    
    % Interpolate Z and color values onto the grid
    Zgrid = griddata(xx, yy, zz, Xgrid, Ygrid, 'natural');
    Cgrid = griddata(xx, yy, cc, Xgrid, Ygrid, 'natural');
    
    % Plot as surface
    surf(Xgrid, Ygrid, Zgrid, Cgrid, 'EdgeColor','none');
    colormap(gca, vmi_cmap);
    colorbar;
    clim([-0.2,0.2])
    zlim([-10000,-4000])
    xlabel('Anterior-Posterior Axis');
    ylabel('Medial-Lateral Axis');
    zlabel('Depth (um)');
    prettyFig;
end

for i = 3:4
    ax1(i) = nexttile;
    monkey = 'scrappy';
    probe = i-2;
    
    cc = CC.SPI(CC.monkey==monkey & CC.probe_index==probe);
    zz = cellfun(@(q) q(2), CC.unit_locations(CC.monkey==monkey & CC.probe_index==probe), 'uni', 1);
    zz = zz - (CC.probe_depth_mm(CC.monkey==monkey & CC.probe_index==probe) * 1000);
    
    xy = CC.probe_gridHole(CC.monkey==monkey & CC.probe_index==probe); 
    
    N = numel(xy);
    out = zeros(N,2);
    for i = 1:N
        vals = xy{i};  % 2x1 cell array of strings
    
        % --- reorder so ML/m/l is first, AP/p/a is second
        if contains(vals{1}, {'AP','p','a'}, 'IgnoreCase', true)
            vals = flipud(vals); % swap order
        end
    
        % --- convert each string to number
        for j = 1:2
            s = vals{j};
    
            if strcmpi(s,'ML') || strcmpi(s,'AP')
                num = 0;
            else
                prefix = lower(s(1));
                val = str2double(s(2:end));
                if prefix == 'm'
                    num = val;
                elseif prefix == 'l'
                    num = -val;
                elseif prefix == 'a'
                    num = val;
                elseif prefix == 'p'
                    num = -val;
                else
                    num = NaN; % fallback if unexpected string
                end
            end
            out(i,j) = num;
        end
    end
    
    xx = out(:,2);
    yy = out(:,1);
    
    % Define grid to interpolate onto
    xlin = linspace(min(xx), max(xx), 50);   % 50 points along X
    ylin = linspace(min(yy), max(yy), 50);   % 50 points along Y
    [Xgrid, Ygrid] = meshgrid(xlin, ylin);
    
    % Interpolate Z and color values onto the grid
    Zgrid = griddata(xx, yy, zz, Xgrid, Ygrid, 'natural');
    Cgrid = griddata(xx, yy, cc, Xgrid, Ygrid, 'natural');
    
    % Plot as surface
    surf(Xgrid, Ygrid, Zgrid, Cgrid, 'EdgeColor','none');
    colormap(gca, spi_cmap);
    colorbar;
    clim([-0.3,0.3])
    zlim([-10000,-4000])
    xlabel('Anterior-Posterior Axis');
    ylabel('Medial-Lateral Axis');
    zlabel('Depth (um)');
    prettyFig;
end



%
X = [xx, yy, zz];

Y = CC.VMI(CC.monkey==monkey & CC.probe_index==probe); 

mdl_vmi = fitlm(X,Y);
mdl_vmi_coeffs = mdl_vmi.Coefficients;
mdl_vmi_anova = anova(mdl_vmi,'summary');

Y = CC.SPI(CC.monkey==monkey & CC.probe_index==probe);

mdl_spi = fitlm(X,Y);
mdl_spi_coeffs = mdl_spi.Coefficients;
mdl_spi_anova = anova(mdl_spi,'summary');






