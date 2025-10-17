%% --------------------------------- KKN PROPOSAL PLOTS --------------------------------- 
clear; clc;
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions'))
addpath(genpath('/Users/kendranoneman/Packages/pairwise_population_metrics'))

%load('/Users/kendranoneman/Data/sapu_dualhemi/concat_tables.mat', 'C', 'Tmdir', 'Tpurs', 'Trfmp');
load('/Users/kendranoneman/Data/sapu_dualhemi/pairs_table.mat');

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

Tmdir.sess_name = categorical(regexprep(string(Tmdir.sess_name),'_g0$',''));
Tpurs.sess_name = categorical(regexprep(string(Tpurs.sess_name),'_g0$',''));
Trfmp.sess_name = categorical(regexprep(string(Trfmp.sess_name),'_g0$',''));
C.sess_name = categorical(regexprep(string(C.sess_name),'_g0$',''));

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

fprintf('Visual: mean = %.3f, std = %.3f, t(%d)=%.2f, p=%.3g\n', mean(a), std(a), statsA.df, statsA.tstat, pA);
fprintf('Motor:  mean = %.3f, std = %.3f, t(%d)=%.2f, p=%.3g\n', mean(b), std(b), statsB.df, statsB.tstat, pB);
fprintf('Pursuit: mean = %.3f, std = %.3f, t(%d)=%.2f, p=%.3g\n', mean(c), std(c), statsC.df, statsC.tstat, pC);

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

%% Distributions of dprime

f3e = figure;
f3e.Position = [100 100 1500 700];
tl = tiledlayout(1,3);
tl.TileSpacing = 'loose';
tl.Padding = 'loose';

ax1(1) = nexttile;
values = CC.dp_vis;
histStyle_KKN(values, 'BIN_WIDTH', 0.25, 'X_LIMITS', [-4,4], 'TEXT_X_POS',2); %, 'Y_LIMITS', [0 600], 'X_LIMITS', [-1,1]);
title('visual epoch')

ax1(2) = nexttile;
values = CC.dp_sac;
histStyle_KKN(values, 'BIN_WIDTH', 0.25, 'X_LIMITS', [-4,4], 'TEXT_X_POS',2); %, 'Y_LIMITS', [0 600], 'X_LIMITS', [-1,1]);
title('saccade epoch')

ax1(3) = nexttile;
values = CC.dp_pur;
histStyle_KKN(values, 'BIN_WIDTH', 0.25, 'X_LIMITS', [-4,4], 'TEXT_X_POS',2); %, 'Y_LIMITS', [0 110], 'X_LIMITS', [-1,1]);
title('pursuit epoch')

linkaxes(ax1,'xy')
xlabel(tl,"d' sensitivity measure (between evoked and spontaneous activity)", 'fontsize', 18)
ylabel(tl,"number of neurons", 'fontsize', 18)

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


%% RSC calculation

nBins = 5;
CC.VMI_quantile = discretize(CC.VMI,  quantile(CC.VMI, linspace(0, 1, nBins+1)));
CC.VMI_quantile(CC.VMI == max(CC.VMI)) = nBins;

CC.SPI_quantile = discretize(CC.SPI, quantile(CC.SPI, linspace(0, 1, nBins+1)));
CC.SPI_quantile(CC.SPI == max(CC.SPI)) = nBins;

sess = unique(CC.sess_name);
pairs_all = cell(length(sess),1);

for s = 18:length(sess)
    fprintf('%d of %d\n', s, numel(sess));
    this_mdir = Tmdir(Tmdir.sess_name==sess(s),:);
    this_purs = Tpurs(Tpurs.sess_name==sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);

    vis_frs_l = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w(1)+50 & u<w(1)+250))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_1, this_mdir.TARG_ON, 'uni', 0);
    vis_frs_r = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w(1)+50 & u<w(1)+250))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_2, this_mdir.TARG_ON, 'uni', 0);

    sac_frs_l = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_1, num2cell(this_mdir.SACCADE), 'uni', 0);
    sac_frs_r = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_2, num2cell(this_mdir.SACCADE), 'uni', 0);

    pur_frs_l = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_purs.spiketimes_1, num2cell(this_purs.pursuitOnset), 'uni', 0);
    pur_frs_r = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_purs.spiketimes_2, num2cell(this_purs.pursuitOnset), 'uni', 0);
    
    pairs_probes = cell(4,1);
    for p1 = 1:2
        if p1==1
            vis_frs_p1 = vis_frs_l;
            sac_frs_p1 = sac_frs_l;
            pur_frs_p1 = pur_frs_l;
        else
            vis_frs_p1 = vis_frs_r;
            sac_frs_p1 = sac_frs_r;
            pur_frs_p1 = pur_frs_r;
        end

        for p2 = 1:2
            if p2==1
                vis_frs_p2 = vis_frs_l;
                sac_frs_p2 = sac_frs_l;
                pur_frs_p2 = pur_frs_l;
            else
                vis_frs_p2 = vis_frs_r;
                sac_frs_p2 = sac_frs_r;
                pur_frs_p2 = pur_frs_r;
            end

            fprintf('%d - %d\n', p1, p2);

            p1_units = CC(CC.sess_name==sess(s) & CC.probe_index==p1,:);
            p2_units = CC(CC.sess_name==sess(s) & CC.probe_index==p2,:);
        
            if p1==p2
                pairs = cell((height(p1_units)*height(p2_units)) - height(p1_units),29);
            else
                pairs = cell(height(p1_units)*height(p2_units),29);
            end

            rr = 1;
            for n1 = 1:height(p1_units)
                for n2 = 1:height(p2_units)  
                    if ~(p1==p2 && n1==n2)
                        n1_clust = p1_units.cluster_id(n1);
                        n2_clust = p2_units.cluster_id(n2);
    
                        % SACCADE TASK
                        TT = this_mdir; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        angs = sort(unique(TT.angle));
                        amps = sort(unique(TT.distance)); 
    
                        %------------------------------ VISUAL ------------------------------%
                        frs_p1 = vis_frs_p1; frs_p2 = vis_frs_p2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if p1_units.vis_sel_dir(n1) > p2_units.vis_sel_dir(n2) %%%%%%%%%%%%%%%%%%
                            bestDir = p1_units.vis_pref_dir(n1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            bestDir = p2_units.vis_pref_dir(n2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
    
                        relAngs = mod(angs - bestDir + 180, 360) - 180;
                        [~, sortIdx] = sort(relAngs);
    
                        [rsc_perDir, fr_perDir] = deal(zeros(1, numel(angs)));
                        n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                        for d = 1:numel(angs)
                            thisAng = angs(d);
    
                            n12_fr = []; n1_zfr = []; n2_zfr = [];
                            for a = 1:numel(amps)
                                n1_fr = cellfun(@(q) q{n1_clust+1}, frs_p1(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
                                n2_fr = cellfun(@(q) q{n2_clust+1}, frs_p2(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
    
                                n12_fr = [n12_fr; n1_fr; n2_fr];
    
                                n1_zfr = [n1_zfr; zscore(n1_fr)];
                                n2_zfr = [n2_zfr; zscore(n2_fr)];
                            end
                            n1_zfr_all = [n1_zfr_all; n1_zfr];
                            n2_zfr_all = [n2_zfr_all; n2_zfr];
                            n12_fr_all = [n12_fr_all; n12_fr];
    
                            [rho,~] = corr(n1_zfr, n2_zfr);
    
                            rsc_perDir(sortIdx(d)) = rho;
                            fr_perDir(sortIdx(d)) = mean(n12_fr);
                        end
    
                        [rsig,~] = corr(cellfun(@(q) q{n1_clust+1}, frs_p1, 'uni', 1), cellfun(@(q) q{n2_clust+1}, frs_p2, 'uni', 1));
                        [rsc,~] = corr(n1_zfr_all, n2_zfr_all);
                        mnFR = mean(n12_fr_all);
    
                        rsig_vis = rsig; rsc_vis = rsc; mnFR_vis = mnFR; rsc_perDir_vis = rsc_perDir; fr_perDir_vis = fr_perDir; 
    
                        %------------------------------ SACCADE ------------------------------%
                        frs_p1 = sac_frs_p1; frs_p2 = sac_frs_p2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if p1_units.sac_sel_dir(n1) > p2_units.sac_sel_dir(n2) %%%%%%%%%%%%%%%%%%
                            bestDir = p1_units.sac_pref_dir(n1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            bestDir = p2_units.sac_pref_dir(n2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
    
                        relAngs = mod(angs - bestDir + 180, 360) - 180;
                        [~, sortIdx] = sort(relAngs);
    
                        [rsc_perDir, fr_perDir] = deal(zeros(1, numel(angs)));
                        n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                        for d = 1:numel(angs)
                            thisAng = angs(d);
    
                            n12_fr = []; n1_zfr = []; n2_zfr = [];
                            for a = 1:numel(amps)
                                n1_fr = cellfun(@(q) q{n1_clust+1}, frs_p1(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
                                n2_fr = cellfun(@(q) q{n2_clust+1}, frs_p2(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
    
                                n12_fr = [n12_fr; n1_fr; n2_fr];
    
                                n1_zfr = [n1_zfr; zscore(n1_fr)];
                                n2_zfr = [n2_zfr; zscore(n2_fr)];
                            end
                            n1_zfr_all = [n1_zfr_all; n1_zfr];
                            n2_zfr_all = [n2_zfr_all; n2_zfr];
                            n12_fr_all = [n12_fr_all; n12_fr];
    
                            [rho,~] = corr(n1_zfr, n2_zfr);
    
                            rsc_perDir(sortIdx(d)) = rho;
                            fr_perDir(sortIdx(d)) = mean(n12_fr);
                        end
    
                        [rsig,~] = corr(cellfun(@(q) q{n1_clust+1}, frs_p1, 'uni', 1), cellfun(@(q) q{n2_clust+1}, frs_p2, 'uni', 1));
                        [rsc,~] = corr(n1_zfr_all, n2_zfr_all);
                        mnFR = mean(n12_fr_all);
    
                        rsig_sac = rsig; rsc_sac = rsc; mnFR_sac = mnFR; rsc_perDir_sac = rsc_perDir; fr_perDir_sac = fr_perDir;
    
                        %------------------------------ PURSUIT ------------------------------%
                        TT = this_purs;
                        angs = sort(unique(TT.angle));
                        spes = sort(unique(TT.pursuitSpeed));
    
                        frs_p1 = pur_frs_p1; frs_p2 = pur_frs_p2; 
                        if p1_units.pur_sel_dir(n1) > p2_units.pur_sel_dir(n2) 
                            bestDir = p1_units.pur_pref_dir(n1);
                        else
                            bestDir = p2_units.pur_pref_dir(n2); 
                        end
    
                        relAngs = mod(angs - bestDir + 180, 360) - 180;
                        [~, sortIdx] = sort(relAngs);
    
                        [rsc_perDir, fr_perDir] = deal(zeros(1, numel(angs)));
                        n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                        for d = 1:numel(angs)
                            thisAng = angs(d);
    
                            n12_fr = []; n1_zfr = []; n2_zfr = [];
                            for a = 1:numel(spes)
                                n1_fr = cellfun(@(q) q{n1_clust+1}, frs_p1(TT.angle==angs(d) & TT.pursuitSpeed==spes(a)), 'uni', 1);
                                n2_fr = cellfun(@(q) q{n2_clust+1}, frs_p2(TT.angle==angs(d) & TT.pursuitSpeed==spes(a)), 'uni', 1);
    
                                n12_fr = [n12_fr; n1_fr; n2_fr];
    
                                n1_zfr = [n1_zfr; zscore(n1_fr)];
                                n2_zfr = [n2_zfr; zscore(n2_fr)];
                            end
                            n1_zfr_all = [n1_zfr_all; n1_zfr];
                            n2_zfr_all = [n2_zfr_all; n2_zfr];
                            n12_fr_all = [n12_fr_all; n12_fr];
    
                            [rho,~] = corr(n1_zfr, n2_zfr);
    
                            rsc_perDir(sortIdx(d)) = rho;
                            fr_perDir(sortIdx(d)) = mean(n12_fr);
                        end
    
                        [rsig,~] = corr(cellfun(@(q) q{n1_clust+1}, frs_p1, 'uni', 1), cellfun(@(q) q{n2_clust+1}, frs_p2, 'uni', 1));
                        [rsc,~] = corr(n1_zfr_all, n2_zfr_all);
                        mnFR = mean(n12_fr_all);
    
                        rsig_pur = rsig; rsc_pur = rsc; mnFR_pur = mnFR; rsc_perDir_pur = rsc_perDir; fr_perDir_pur = fr_perDir;
    
                        %---------------------------------------------------------------------------------------------------------
                        
                        pairs(rr,:) = {p1_units.monkey(n1), p1_units.sess_name(n1) ...
                                p1, p2, n1_clust, n2_clust, ...
                                p1_units.VMI_quantile(n1), p2_units.VMI_quantile(n2), ...
                                p1_units.SPI_quantile(n1), p2_units.SPI_quantile(n2), ...
                                norm(p1_units.unit_locations{n1}(1:2) - p2_units.unit_locations{n2}(1:2)), ...
                                abs(mod(p1_units.vis_pref_dir(n1) - p2_units.vis_pref_dir(n2) + 180, 360) - 180), ...
                                abs(mod(p1_units.sac_pref_dir(n1) - p2_units.sac_pref_dir(n2) + 180, 360) - 180), ...
                                abs(mod(p1_units.pur_pref_dir(n1) - p2_units.pur_pref_dir(n2) + 180, 360) - 180), ...
                                mnFR_vis, rsc_vis, rsig_vis, rsc_perDir_vis, fr_perDir_vis, ...
                                mnFR_sac, rsc_sac, rsig_sac, rsc_perDir_sac, fr_perDir_sac, ...
                                mnFR_pur, rsc_pur, rsig_pur, rsc_perDir_pur, fr_perDir_pur
                                };
                        rr = rr + 1;
                    end
                end
            end   
            pairs_probes{(p1-1)*2 + p2} = pairs;
        end
    end
    pairs_all{s} = pairs_probes;
end

pairs_all2 = vertcat(pairs_all{:});
pairs_all3 = vertcat(pairs_all2{:});
pairs_tbl = cell2table(pairs_all3,'VariableNames',{'monkey','sess_name','n1_probe','n2_probe','n1_cluster','n2_cluster','n1_VMI','n2_VMI','n1_SPI','n2_SPI','depth_diff','vis_prefDir_diff','sac_prefDir_diff','pur_prefDir_diff','mnFR_vis','rsc_vis','rsig_vis','rsc_perDir_vis','fr_perDir_vis','mnFR_sac','rsc_sac','rsig_sac','rsc_perDir_sac','fr_perDir_sac','mnFR_pur','rsc_pur','rsig_pur','rsc_perDir_pur','fr_perDir_pur'});

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

%% CALCULATE POP MATRIX

sess = unique(CC.sess_name);
zDims = 0:10;
num_repeats = 10;

all_pop_mats = nan(length(sess),3,3,5);
for s = 1:length(sess)
    if s~=5
        fprintf('%d of %d\n', s, numel(sess));
        this_mdir = Tmdir(Tmdir.sess_name==sess(s),:);
        this_purs = Tpurs(Tpurs.sess_name==sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
    
        frs_l = cellfun(@(q,w) cellfun(@(u) (sum(u>=w(1)+50 & u<w(1)+250)), q, 'uni', 1), this_mdir.spiketimes_1, this_mdir.TARG_ON, 'uni', 0);
        frs_r = cellfun(@(q,w) cellfun(@(u) (sum(u>=w(1)+50 & u<w(1)+250)), q, 'uni', 1), this_mdir.spiketimes_2, this_mdir.TARG_ON, 'uni', 0);
    
        frs_l = vertcat(frs_l{:}); frs_r = vertcat(frs_r{:});
        frs_l_vis = frs_l(:,CC.cluster_id(CC.sess_name==sess(s) & CC.probe_index==1)+1);
        frs_r_vis = frs_r(:,CC.cluster_id(CC.sess_name==sess(s) & CC.probe_index==2)+1);
    
        frs_l = cellfun(@(q,w) cellfun(@(u) (sum(u>=w-100 & u<w+100)), q, 'uni', 1), this_mdir.spiketimes_1, num2cell(this_mdir.SACCADE), 'uni', 0);
        frs_r = cellfun(@(q,w) cellfun(@(u) (sum(u>=w-100 & u<w+100)), q, 'uni', 1), this_mdir.spiketimes_2, num2cell(this_mdir.SACCADE), 'uni', 0);
    
        frs_l = vertcat(frs_l{:}); frs_r = vertcat(frs_r{:});
        frs_l_sac = frs_l(:,CC.cluster_id(CC.sess_name==sess(s) & CC.probe_index==1)+1);
        frs_r_sac = frs_r(:,CC.cluster_id(CC.sess_name==sess(s) & CC.probe_index==2)+1);
    
        frs_l = cellfun(@(q,w) cellfun(@(u) (sum(u>=w-100 & u<w+100)), q, 'uni', 1), this_purs.spiketimes_1, num2cell(this_purs.pursuitOnset), 'uni', 0);
        frs_r = cellfun(@(q,w) cellfun(@(u) (sum(u>=w-100 & u<w+100)), q, 'uni', 1), this_purs.spiketimes_2, num2cell(this_purs.pursuitOnset), 'uni', 0);
        
        frs_l = vertcat(frs_l{:}); frs_r = vertcat(frs_r{:});
        frs_l_pur = frs_l(:,CC.cluster_id(CC.sess_name==sess(s) & CC.probe_index==1)+1);
        frs_r_pur = frs_r(:,CC.cluster_id(CC.sess_name==sess(s) & CC.probe_index==2)+1);
    
        [r_idx,l_idx] = deal(cell(1,num_repeats));
        for i = 1:num_repeats
            r_idx{i} = sort(randperm(size(frs_r_vis,2),(floor(min(size(frs_l_vis,2),size(frs_r_vis,2))/2))));
            l_idx{i} = sort(randperm(size(frs_l_vis,2),floor(min(size(frs_l_vis,2),size(frs_r_vis,2))/2)));
        end
    
        % SACCADE TASK
        TT = this_mdir; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        angs = sort(unique(TT.angle));
    
        % VIS
        frs_l = frs_l_vis;
        frs_r = frs_r_vis;
    
        X_l = []; X_r = [];
        for d = 1:numel(angs)
            X_l = [X_l; zscore(frs_l(TT.angle==angs(d),:))];
            X_r = [X_r; zscore(frs_l(TT.angle==angs(d),:))];
        end
    
        [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X_l');
        [~, ls, ~, psv, dim, ~] = compute_population_metrics(X_l', zDims);
        if ~isempty(ls)
            all_pop_mats(s,1,1,:) = [rsc_mean, rsc_std, ls(1), psv, dim];
        end
    
        [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X_r');
        [~, ls, ~, psv, dim, ~] = compute_population_metrics(X_r', zDims);
        if ~isempty(ls)
            all_pop_mats(s,1,2,:) = [rsc_mean, rsc_std, ls(1), psv, dim];
        end
    
        [rsc_mean_all, rsc_std_all, ls_all, psv_all, dim_all] = deal(nan(1,num_repeats));
        for i = 1:num_repeats
            X = []; 
            for d = 1:numel(angs)
                X = [X; zscore([frs_l(TT.angle==angs(d),l_idx{i}) frs_r(TT.angle==angs(d),r_idx{i})])];
            end
    
            [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X');
            [~, ls, ~, psv, dim, ~] = compute_population_metrics(X', zDims);
    
            rsc_mean_all(i) = rsc_mean; rsc_std_all(i) = rsc_std;
            if ~isempty(ls)
                ls_all(i) = ls(1); psv_all(i) = psv; dim_all(i) = dim;
            end
        end
    
        all_pop_mats(s,1,3,:) = [mean(rsc_mean_all,'omitnan'), mean(rsc_std_all,'omitnan'), mean(ls_all,'omitnan'), mean(psv_all,'omitnan'), mean(dim_all,'omitnan')];
    
        % MDIR
        frs_l = frs_l_sac;
        frs_r = frs_r_sac;
    
        X_l = []; X_r = [];
        for d = 1:numel(angs)
            X_l = [X_l; zscore(frs_l(TT.angle==angs(d),:))];
            X_r = [X_r; zscore(frs_l(TT.angle==angs(d),:))];
        end
    
        [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X_l');
        [~, ls, ~, psv, dim, ~] = compute_population_metrics(X_l', zDims);
        if ~isempty(ls)
            all_pop_mats(s,2,1,:) = [rsc_mean, rsc_std, ls(1), psv, dim];
        end
    
        [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X_r');
        [~, ls, ~, psv, dim, ~] = compute_population_metrics(X_r', zDims);
        if ~isempty(ls)
            all_pop_mats(s,2,2,:) = [rsc_mean, rsc_std, ls(1), psv, dim];
        end
    
        [rsc_mean_all, rsc_std_all, ls_all, psv_all, dim_all] = deal(nan(1,num_repeats));
        for i = 1:num_repeats
            X = []; 
            for d = 1:numel(angs)
                X = [X; zscore([frs_l(TT.angle==angs(d),l_idx{i}) frs_r(TT.angle==angs(d),r_idx{i})])];
            end
    
            [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X');
            [~, ls, ~, psv, dim, ~] = compute_population_metrics(X', zDims);
    
            rsc_mean_all(i) = rsc_mean; rsc_std_all(i) = rsc_std;
            if ~isempty(ls)
                ls_all(i) = ls(1); psv_all(i) = psv; dim_all(i) = dim;
            end
        end
    
        all_pop_mats(s,2,3,:) = [mean(rsc_mean_all,'omitnan'), mean(rsc_std_all,'omitnan'), mean(ls_all,'omitnan'), mean(psv_all,'omitnan'), mean(dim_all,'omitnan')];
    
    
        % PURSUIT TASK
        TT = this_purs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        angs = sort(unique(TT.angle));
        spes = sort(unique(TT.pursuitSpeed)); 
    
        frs_l = frs_l_pur;
        frs_r = frs_r_pur;
    
        X_l = []; X_r = [];
        for d = 1:numel(angs)
            X_l = [X_l; zscore(frs_l(TT.angle==angs(d),:))];
            X_r = [X_r; zscore(frs_l(TT.angle==angs(d),:))];
        end
    
        [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X_l');
        [~, ls, ~, psv, dim, ~] = compute_population_metrics(X_l', zDims);
        if ~isempty(ls)
            all_pop_mats(s,3,1,:) = [rsc_mean, rsc_std, ls(1), psv, dim];
        end
    
        [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X_r');
        [~, ls, ~, psv, dim, ~] = compute_population_metrics(X_r', zDims);
        if ~isempty(ls)
            all_pop_mats(s,3,2,:) = [rsc_mean, rsc_std, ls(1), psv, dim];
        end
    
        [rsc_mean_all, rsc_std_all, ls_all, psv_all, dim_all] = deal(nan(1,num_repeats));
        for i = 1:num_repeats
            X = []; 
            for d = 1:numel(angs)
                X = [X; zscore([frs_l(TT.angle==angs(d),l_idx{i}) frs_r(TT.angle==angs(d),r_idx{i})])];
            end
    
            [~, rsc_mean, rsc_std] = compute_pairwise_metrics(X');
            [~, ls, ~, psv, dim, ~] = compute_population_metrics(X', zDims);
    
            rsc_mean_all(i) = rsc_mean; rsc_std_all(i) = rsc_std;
            if ~isempty(ls)
                ls_all(i) = ls(1); psv_all(i) = psv; dim_all(i) = dim;
            end
        end
    
        all_pop_mats(s,3,3,:) = [mean(rsc_mean_all,'omitnan'), mean(rsc_std_all,'omitnan'), mean(ls_all,'omitnan'), mean(psv_all,'omitnan'), mean(dim_all,'omitnan')];
    end
end

%% Remove pairs of same neurons
pairs_tbl(pairs_tbl.n1_probe==pairs_tbl.n2_probe & pairs_tbl.n1_cluster==pairs_tbl.n2_cluster, :) = [];
pairs_tbl(pairs_tbl.rsc_vis==1 | pairs_tbl.rsc_sac==1 | pairs_tbl.rsc_pur==1,:) = [];
pairs_tbl(pairs_tbl.depth_diff==0,:) = [];
pairs_tbl(pairs_tbl.vis_prefDir_diff==0 | pairs_tbl.sac_prefDir_diff==0 | pairs_tbl.pur_prefDir_diff==0,:) = [];
pairs_tbl_uni = pairs_tbl((pairs_tbl.n1_cluster < pairs_tbl.n2_cluster) | (pairs_tbl.n1_probe < pairs_tbl.n2_probe),:);

% Count up number of pairs within each hemi and across hemi
within_hemi_pairs = sum(pairs_tbl_uni.n1_probe==pairs_tbl_uni.n2_probe);
across_hemi_pairs = sum(pairs_tbl_uni.n1_probe~=pairs_tbl_uni.n2_probe); 

%% RSC as a function of distance and preferred direction diff (within hemi)
P = pairs_tbl_uni(pairs_tbl_uni.n1_probe==pairs_tbl_uni.n2_probe,:);

f3f = figure;
f3f.Position = [100 100 800 400];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

N_BINS = 15;

% DISTANCE
ax1(1) = nexttile;

% even width 
%bin_edges = linspace(0,7000,N_BINS+1);
bin_edges = 0:500:7500;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

[mean_rho,sem_rho,n_pairs] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_vis(P.depth_diff >= bin_edges(b) & P.depth_diff < bin_edges(b+1)));
    n_pairs(b) = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs(b));
end
errorbar(bin_centers./1000, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', green_visual./255);
hold on;

[r, pval] = corr(P.depth_diff, rtoZ(P.rsc_vis), 'Type', 'Spearman');

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_sac(P.depth_diff >= bin_edges(b) & P.depth_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers./1000, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', purple_saccadey./255);
hold on;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_pur(P.depth_diff >= bin_edges(b) & P.depth_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers./1000, mean_rho, sem_rho, 'go:', 'LineWidth', 2, 'Color', red_pursuity./255);

prettyFig;
axis square

% PREF DIR
ax1(2) = nexttile;

% even width 
%bin_edges = linspace(,7000,N_BINS+1);
bin_edges = 0:20:180;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_vis(P.vis_prefDir_diff >= bin_edges(b) & P.vis_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', green_visual./255);
hold on;

[r, pval] = corr(P.vis_prefDir_diff, rtoZ(P.rsc_vis), 'Type', 'Spearman');

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_sac(P.sac_prefDir_diff>= bin_edges(b) & P.sac_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', purple_saccadey./255);
hold on;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_pur(P.pur_prefDir_diff >= bin_edges(b) & P.pur_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'go:', 'LineWidth', 2, 'Color', red_pursuity./255);

xticks(bin_centers)
xlim([5,175])

prettyFig;
axis square

linkaxes(ax1,'y')

savebigPDF(f3f, '/Users/kendranoneman/Milestones/proposal/figs/rsc_tuning.pdf') 


%% Heatmaps of rsc, VMI quantiles (within hemi)
P = pairs_tbl(pairs_tbl.n1_probe==pairs_tbl.n2_probe,:);

f3g = figure;
f3g.Position = [100 100 800 400]; 
tl = tiledlayout(1,3); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'compact'; 

xx = P.n1_VMI;
yy = P.n2_VMI;

% VISUAL
ax1(1) = nexttile;
cc = P.rsc_vis;

% Compute mean values per (xx, yy)
M1 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);

% Compute counts per (xx, yy)
N1 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M1);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(green_visual,0,256);
clim([0.0016,0.0282]);
colormap(ax1(1), cmap);
cb = colorbar(ax1(1));

prettyFig;

% SACCADE
ax1(2) = nexttile;
cc = P.rsc_sac;

M2 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N2 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M2);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(purple_saccadey,0,256);
clim([0.0016,0.0282]);
colormap(ax1(2), cmap);
cb = colorbar(ax1(2));

prettyFig;

% VISUAL
ax1(3) = nexttile;
cc = P.rsc_pur;

M3 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N3 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M3);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(red_pursuity,0,256);
clim([0.0016,0.0282]);
colormap(ax1(3), cmap);
cb = colorbar(ax1(3));

prettyFig;

%savebigPDF(f3g, '/Users/kendranoneman/Milestones/proposal/figs/rsc_heatmap_vmi.pdf') 

%% Heatmaps of rsc, SPI quantiles (within hemi)
P = pairs_tbl(pairs_tbl.n1_probe==pairs_tbl.n2_probe,:);

f3h = figure;
f3h.Position = [100 100 800 400]; 
tl = tiledlayout(1,3); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'compact'; 

xx = P.n1_SPI;
yy = P.n2_SPI;

% VISUAL
ax1(1) = nexttile;
cc = P.rsc_vis;

% Compute mean values per (xx, yy)
M4 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);

% Compute counts per (xx, yy)
N4 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M4);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(green_visual,0,256);
clim([0.0016,0.0282]);
colormap(ax1(1), cmap);
cb = colorbar(ax1(1));

prettyFig;

% SACCADE
ax1(2) = nexttile;
cc = P.rsc_sac;

M5 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N5 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M5);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(purple_saccadey,0,256);
clim([0.0016,0.0282]);
colormap(ax1(2), cmap);
cb = colorbar(ax1(2));

prettyFig;

% VISUAL
ax1(3) = nexttile;
cc = P.rsc_pur;

M6 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N6 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M6);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(red_pursuity,0,256);
clim([0.0016,0.0282]);
colormap(ax1(3), cmap);
cb = colorbar(ax1(3));

prettyFig;
savebigPDF(f3h, '/Users/kendranoneman/Milestones/proposal/figs/rsc_heatmap_spi.pdf') 




%% CHAPTER 4 - SACCADE-PURSUIT INTERACTIONS ACROSS HEMIS

%% PSTH colored by VMI and SPI quantiles (contra, ipsi)
nBins = 5;
sigma = 10;
visWindow = [-100,300];
sacWindow = [-200,200];
purWindow = [-200,200];

CC.VMI_quantile = discretize(CC.VMI,  quantile(CC.VMI, linspace(0, 1, nBins+1)));
CC.VMI_quantile(CC.VMI == max(CC.VMI)) = nBins;

[vis_mn_vmi1,vis_sem_vmi1,sac_mn_vmi1,sac_sem_vmi1,pur_mn_vmi1,pur_sem_vmi1] = deal(cell(numel(these_sess),nBins));
[vis_mn_vmi2,vis_sem_vmi2,sac_mn_vmi2,sac_sem_vmi2,pur_mn_vmi2,pur_sem_vmi2] = deal(cell(numel(these_sess),nBins));
for q = 1:nBins
    fprintf('Bin %d of %d...\n', q, nBins);
    these_units = CC(CC.VMI_quantile==q,:);

    these_sess = unique(these_units.sess_name);
    for s = 1:numel(these_sess)
        this_mdirTbl = Tmdir(Tmdir.sess_name==these_sess(s),:);
        this_pursTbl = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
        
        these_spikes1a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_mdirTbl.spiketimes_1, 'uni', 0);
        these_spikes2a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_mdirTbl.spiketimes_2, 'uni', 0);

        for i = 1:2
            % VISUAL PSTH
            time_window = visWindow;

            spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes1a, this_mdirTbl.TARG_ON, 'uni', 0);
            spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes2a, this_mdirTbl.TARG_ON, 'uni', 0);

            if i == 1 % CONTRA
                spike_times1 = spike_times1(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
                spike_times2 = spike_times2(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
            else
                spike_times1 = spike_times1(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
                these_spikes2 = spike_times2(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
            end
    
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
    
            if i == 1 
                vis_mn_vmi1{s,q} = vis_mn;
                vis_sem_vmi1{s,q} = vis_sem;
            else
                vis_mn_vmi2{s,q} = vis_mn;
                vis_sem_vmi2{s,q} = vis_sem;
            end
        end

        for i = 1:2
            % SAC PSTH
            time_window = sacWindow;
    
            spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1a, num2cell(this_mdirTbl.SACCADE), 'uni', 0);
            spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2a, num2cell(this_mdirTbl.SACCADE), 'uni', 0);
    
            if i == 1 % CONTRA
                spike_times1 = spike_times1(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
                spike_times2 = spike_times2(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
            else
                spike_times1 = spike_times1(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
                these_spikes2 = spike_times2(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
            end
    
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

            if i == 1 
                sac_mn_vmi1{s,q} = sac_mn;
                sac_sem_vmi1{s,q} = sac_sem;
            else
                sac_mn_vmi2{s,q} = sac_mn;
                sac_sem_vmi2{s,q} = sac_sem;
            end
        end

        % PURS PSTH
        time_window = purWindow;

        these_spikes1a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl.spiketimes_1, 'uni', 0);
        these_spikes2a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl.spiketimes_2, 'uni', 0);

        for i = 1:2

            spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1a, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);
            spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2a, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);
    
            if i == 1 % CONTRA
                spike_times1 = spike_times1(ismember(this_pursTbl.angle,[0,45,90,270,315]));
                spike_times2 = spike_times2(ismember(this_pursTbl.angle,[190,135,180,225,315]));
            else
                spike_times1 = spike_times1(ismember(this_pursTbl.angle,[190,135,180,225,315]));
                these_spikes2 = spike_times2(ismember(this_pursTbl.angle,[0,45,90,270,315]));
            end
    
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
    
            if i == 1
                pur_mn_vmi1{s,q} = pur_mn;
                pur_sem_vmi1{s,q} = pur_sem;
            else
                pur_mn_vmi2{s,q} = pur_mn;
                pur_sem_vmi2{s,q} = pur_sem;
            end
        end
    end
end

CC.SPI_quantile = discretize(CC.SPI, quantile(CC.SPI, linspace(0, 1, nBins+1)));
CC.SPI_quantile(CC.SPI == max(CC.SPI)) = nBins;

[vis_mn_spi1,vis_sem_spi1,sac_mn_spi1,sac_sem_spi1,pur_mn_spi1,pur_sem_spi1] = deal(cell(numel(these_sess),nBins));
[vis_mn_spi2,vis_sem_spi2,sac_mn_spi2,sac_sem_spi2,pur_mn_spi2,pur_sem_spi2] = deal(cell(numel(these_sess),nBins));
for q = 1:nBins
    fprintf('Bin %d of %d...\n', q, nBins);
    these_units = CC(CC.SPI_quantile==q,:);

    these_sess = unique(these_units.sess_name);
    for s = 1:numel(these_sess)
        this_mdirTbl = Tmdir(Tmdir.sess_name==these_sess(s),:);
        this_pursTbl = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
        
        these_spikes1a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_mdirTbl.spiketimes_1, 'uni', 0);
        these_spikes2a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_mdirTbl.spiketimes_2, 'uni', 0);

        for i = 1:2
            % VISUAL PSTH
            time_window = visWindow;

            spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes1a, this_mdirTbl.TARG_ON, 'uni', 0);
            spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v(1), w, 'uni', 0), these_spikes2a, this_mdirTbl.TARG_ON, 'uni', 0);

            if i == 1 % CONTRA
                these_spikes1 = these_spikes1(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
                these_spikes2 = these_spikes2(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
            else
                these_spikes1 = these_spikes1(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
                these_spikes2 = these_spikes2(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
            end
    
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
    
            if i == 1 
                vis_mn_spi1{s,q} = vis_mn;
                vis_sem_spi1{s,q} = vis_sem;
            else
                vis_mn_spi2{s,q} = vis_mn;
                vis_sem_spi2{s,q} = vis_sem;
            end
        end

        % SAC PSTH
        time_window = sacWindow;

        for i = 1:2
            % SAC PSTH
            time_window = sacWindow;
    
            spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1a, num2cell(this_mdirTbl.SACCADE), 'uni', 0);
            spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2a, num2cell(this_mdirTbl.SACCADE), 'uni', 0);
    
            if i == 1 % CONTRA
                these_spikes1 = these_spikes1(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
                these_spikes2 = these_spikes2(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
            else
                these_spikes1 = these_spikes1(ismember(this_mdirTbl.angle,[190,135,180,225,315]));
                these_spikes2 = these_spikes2(ismember(this_mdirTbl.angle,[0,45,90,270,315]));
            end
    
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

            if i == 1 
                sac_mn_spi1{s,q} = sac_mn;
                sac_sem_spi1{s,q} = sac_sem;
            else
                sac_mn_spi2{s,q} = sac_mn;
                sac_sem_spi2{s,q} = sac_sem;
            end
        end

        % PURS PSTH
        time_window = purWindow;

        these_spikes1a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl.spiketimes_1, 'uni', 0);
        these_spikes2a = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl.spiketimes_2, 'uni', 0);

        for i = 1:2
            spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1a, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);
            spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2a, num2cell(this_pursTbl.pursuitOnset), 'uni', 0);
    
            if i == 1 % CONTRA
                these_spikes1 = these_spikes1(ismember(this_pursTbl.angle,[0,45,90,270,315]));
                these_spikes2 = these_spikes2(ismember(this_pursTbl.angle,[190,135,180,225,315]));
            else
                these_spikes1 = these_spikes1(ismember(this_pursTbl.angle,[190,135,180,225,315]));
                these_spikes2 = these_spikes2(ismember(this_pursTbl.angle,[0,45,90,270,315]));
            end
    
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

            if i == 1
                pur_mn_spi1{s,q} = pur_mn;
                pur_sem_spi1{s,q} = pur_sem;
            else
                pur_mn_spi2{s,q} = pur_mn;
                pur_sem_spi2{s,q} = pur_sem;
            end
        end
    end
end

%% Distribution of preferred directions
bw = 15;
face_alpha = 0.1;

f4a = figure;
f4a.Position = [100 100 600 600];  

% Scrappy 
RLIM = 600;

values = [CC.vis_pref_dir(CC.probe_index==1); wrapTo360(CC.vis_pref_dir(CC.probe_index==2)+180)];
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', green_visual./255, 'RLIM', RLIM, ...
                   'LINE_COLOR', green_visual./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = [CC.sac_pref_dir(CC.probe_index==1); wrapTo360(CC.sac_pref_dir(CC.probe_index==2)+180)];
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', purple_saccadey./255, ...
                  'LINE_COLOR', purple_saccadey./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
hold on;
values = [CC.pur_pref_dir(CC.probe_index==1); wrapTo360(CC.pur_pref_dir(CC.probe_index==2)+180)];
polarHistStyle_KKN(values, 'BIN_WIDTH', bw, 'FACE_COLOR', red_pursuity./255, ...
                  'LINE_COLOR', red_pursuity./255, 'FACE_ALPHA', face_alpha, 'MARK_RAD', RLIM);
prettyFig; 

savebigPDF(f4a, '/Users/kendranoneman/Milestones/proposal/figs/polarhists_prefdir.pdf') 

%% Difference between preferred directions

f4b = figure;
f4b.Position = [100 100 900 600];
tl = tiledlayout(3,1);
tl.TileSpacing = 'tight';
tl.Padding = 'tight';

vis_sac_diff = abs(mod(CC.vis_pref_dir - CC.sac_pref_dir + 180, 360) - 180);
vis_pur_diff = abs(mod(CC.vis_pref_dir - CC.pur_pref_dir + 180, 360) - 180);
sac_pur_diff = abs(mod(CC.sac_pref_dir - CC.pur_pref_dir + 180, 360) - 180);

ax1(1) = nexttile(1);
histStyle_KKN(vis_sac_diff, 'BIN_WIDTH', 10, 'X_LIMITS', [0,180], 'Y_LIMITS', [0 600]);
xticks([0,40,80,120,160])

ax1(2) = nexttile(2);
histStyle_KKN(vis_pur_diff, 'BIN_WIDTH', 10, 'X_LIMITS', [0,180], 'Y_LIMITS', [0 600]);
xticks([0,40,80,120,160])

ax1(3) = nexttile(3);
histStyle_KKN(sac_pur_diff, 'BIN_WIDTH', 10, 'X_LIMITS', [0,180], 'Y_LIMITS', [0 600]);
xticks([0,40,80,120,160])

prettyFig;

savebigPDF(f4b, '/Users/kendranoneman/Milestones/proposal/figs/dists_prefdirdiff.pdf') 

%% Difference between preferred directions as a function of VMI or SPI
nBins = 10;

spi_cmap = customColormap(red_pursuity,purple_saccadey, nBins);
vmi_cmap = customColormap(purple_saccadey, green_visual, nBins);

vis_sac_diff = abs(mod(CC.vis_pref_dir - CC.sac_pref_dir + 180, 360) - 180);
vis_pur_diff = abs(mod(CC.vis_pref_dir - CC.pur_pref_dir + 180, 360) - 180);
sac_pur_diff = abs(mod(CC.sac_pref_dir - CC.pur_pref_dir + 180, 360) - 180);

f4c = figure;
f4c.Position = [100 100 1400 600];
tl = tiledlayout(3,1);

ax1(1) = nexttile;

xx_edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
bin_centers = xx_edges(1:end-1) + diff(xx_edges)/2;
xx = discretize(CC.VMI,  xx_edges);

[mns,sems] = deal(nan(1,nBins));
for x = 1:nBins
    [mn, sem, yu, yl] = sem_errorbar(vis_sac_diff(xx==x));
    mns(x) = mn;
    sems(x) = sem;

    errorbar(x, mn, sem, 'o-', 'LineWidth', 3, 'Color', vmi_cmap(x,:))
end
hold on;
plot(1:numel(bin_centers), mns, 'k-', 'LineWidth', 2);

xx_edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
bin_centers = xx_edges(1:end-1) + diff(xx_edges)/2;
xx = discretize(CC.SPI,  xx_edges);

[mns,sems] = deal(nan(1,nBins));
for x = 1:nBins
    [mn, sem, yu, yl] = sem_errorbar(vis_sac_diff(xx==x));
    mns(x) = mn;
    sems(x) = sem;

    errorbar(x, mn, sem, 'o-', 'LineWidth', 3, 'Color', spi_cmap(x,:))
end
plot(1:numel(bin_centers), mns, 'k-', 'LineWidth', 2);
prettyFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1(2) = nexttile;

xx_edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
bin_centers = xx_edges(1:end-1) + diff(xx_edges)/2;
yy = discretize(CC.VMI,  xx_edges);
xx = discretize(vis_pur_diff, 0:20:180);

[mns,sems] = deal(nan(1,nBins));
for x = 1:nBins
    [mn, sem, yu, yl] = sem_errorbar(vis_pur_diff(xx==x));
    mns(x) = mn;
    sems(x) = sem;
end

errorbar(bin_centers, mns, sems, 'o-', 'LineWidth', 2, 'Color', green_visual./255);
xticks(bin_centers)

hold on;
xx_edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
bin_centers = xx_edges(1:end-1) + diff(xx_edges)/2;
xx = discretize(CC.SPI,  xx_edges);

[mns,sems] = deal(nan(1,nBins));
for x = 1:nBins
    [mn, sem, yu, yl] = sem_errorbar(vis_pur_diff(xx==x));
    mns(x) = mn;
    sems(x) = sem;
end

errorbar(bin_centers, mns, sems, 'o-', 'LineWidth', 2, 'Color', red_pursuity./255);
xticks(bin_centers)
prettyFig;

ax1(3) = nexttile;

xx_edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
bin_centers = xx_edges(1:end-1) + diff(xx_edges)/2;
yy = discretize(CC.VMI,  xx_edges);
xx = discretize(sac_pur_diff, 0:20:180);

[mns,sems] = deal(nan(1,nBins));
for x = 1:nBins
    [mn, sem, yu, yl] = sem_errorbar(sac_pur_diff(xx==x));
    mns(x) = mn;
    sems(x) = sem;
end

errorbar(1:numel(bin_centers), mns, sems, 'o-', 'LineWidth', 2, 'Color', green_visual./255);
xticks(bin_centers)

hold on;
xx_edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
bin_centers = xx_edges(1:end-1) + diff(xx_edges)/2;
xx = discretize(CC.SPI,  xx_edges);

[mns,sems] = deal(nan(1,nBins));
for x = 1:nBins
    [mn, sem, yu, yl] = sem_errorbar(sac_pur_diff(xx==x));
    mns(x) = mn;
    sems(x) = sem;
end

errorbar(1:numel(bin_centers), mns, sems, 'o-', 'LineWidth', 2, 'Color', red_pursuity./255);
xticks(bin_centers)
prettyFig;


%%

f4c = figure;
f4c.Position = [100 100 800 1000];
tl = tiledlayout(3,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

vis_sac_diff = abs(mod(CC.vis_pref_dir - CC.sac_pref_dir + 180, 360) - 180);
vis_pur_diff = abs(mod(CC.vis_pref_dir - CC.pur_pref_dir + 180, 360) - 180);
sac_pur_diff = abs(mod(CC.sac_pref_dir - CC.pur_pref_dir + 180, 360) - 180);

y_edges = 0:20:180;

% vis versus sac 
ax(1) = nexttile(1);
xx_edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
xx = discretize(CC.VMI,  xx_edges);
yy = discretize(vis_sac_diff, y_edges);

counts = histcounts2(xx, yy, min(xx):max(xx)+1, min(yy):max(yy)+1);

imagesc(counts');        
axis xy;                        
colormap('gray'); 
colorbar;
clim([30,140]);

xticks((1:length(xx_edges))-0.5);
xticklabels(compose('%.2f', xx_edges));
xtickangle(0);

y_centers = y_edges(1:end-1) + diff(y_edges)/2;
yticks((1:length(x_centers))+0.5);
yticklabels(compose('%.0f', y_centers+10));
prettyFig;


ax(2) = nexttile(2);
xx_edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
xx = discretize(CC.SPI,  xx_edges);
yy = discretize(vis_sac_diff, y_edges);

counts = histcounts2(xx, yy, min(xx):max(xx)+1, min(yy):max(yy)+1);

imagesc(counts');        
axis xy;                        
colormap('gray'); 
colorbar;
clim([30,140]);

xticks((1:length(xx_edges))-0.5);
xticklabels(compose('%.2f', xx_edges));
xtickangle(0);

y_centers = y_edges(1:end-1) + diff(y_edges)/2;
yticks((1:length(x_centers))+0.5);
yticklabels(compose('%.0f', y_centers+10));
prettyFig;


% vis versus pursuit
ax(3) = nexttile(3);
xx_edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
xx = discretize(CC.VMI,  xx_edges);
yy = discretize(vis_pur_diff, y_edges);

counts = histcounts2(xx, yy, min(xx):max(xx)+1, min(yy):max(yy)+1);

imagesc(counts');        
axis xy;                        
colormap('gray'); 
colorbar;
clim([30,140]);

xticks((1:length(xx_edges))-0.5);
xticklabels(compose('%.2f', xx_edges));
xtickangle(0);

y_centers = y_edges(1:end-1) + diff(y_edges)/2;
yticks((1:length(x_centers))+0.5);
yticklabels(compose('%.0f', y_centers+10));
prettyFig;



ax(4) = nexttile(4);
xx_edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
xx = discretize(CC.SPI,  xx_edges);
yy = discretize(vis_pur_diff, y_edges);

counts = histcounts2(xx, yy, min(xx):max(xx)+1, min(yy):max(yy)+1);

imagesc(counts');        
axis xy;                        
colormap('gray'); 
colorbar;
clim([30,140]);

xticks((1:length(xx_edges))-0.5);
xticklabels(compose('%.2f', xx_edges));
xtickangle(0);

y_centers = y_edges(1:end-1) + diff(y_edges)/2;
yticks((1:length(x_centers))+0.5);
yticklabels(compose('%.0f', y_centers+10));
prettyFig;



% sac versus pursuit
ax(5) = nexttile(5);
xx_edges = quantile(CC.VMI, linspace(0, 1, nBins+1));
xx = discretize(CC.VMI,  xx_edges);
yy = discretize(sac_pur_diff, y_edges);

counts = histcounts2(xx, yy, min(xx):max(xx)+1, min(yy):max(yy)+1);

imagesc(counts');        
axis xy;                        
colormap('gray'); 
colorbar;
clim([30,140]);

xticks((1:length(xx_edges))-0.5);
xticklabels(compose('%.2f', xx_edges));
xtickangle(0);

y_centers = y_edges(1:end-1) + diff(y_edges)/2;
yticks((1:length(x_centers))+0.5);
yticklabels(compose('%.0f', y_centers+10));
prettyFig;


ax(6) = nexttile(6);
xx_edges = quantile(CC.SPI, linspace(0, 1, nBins+1));
xx = discretize(CC.SPI,  xx_edges);
yy = discretize(sac_pur_diff, y_edges);

counts = histcounts2(xx, yy, min(xx):max(xx)+1, min(yy):max(yy)+1);

imagesc(counts');        
axis xy;                        
colormap('gray'); 
colorbar;
clim([30,140]);

xticks((1:length(xx_edges))-0.5);
xticklabels(compose('%.2f', xx_edges));
xtickangle(0);

y_centers = y_edges(1:end-1) + diff(y_edges)/2;
yticks((1:length(x_centers))+0.5);
yticklabels(compose('%.0f', y_centers+10));
prettyFig;

%linkaxes(ax2,'xy')


%% RSC as a function of preferred direction diff (within versus across hemi)
P = pairs_tbl_uni(pairs_tbl_uni.n1_probe==pairs_tbl_uni.n2_probe,:);

f4f = figure;
f4f.Position = [100 100 800 400];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

N_BINS = 15;

% DISTANCE
ax1(1) = nexttile;

% even width 
%bin_edges = linspace(,7000,N_BINS+1);
bin_edges = 0:20:180;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_vis(P.vis_prefDir_diff >= bin_edges(b) & P.vis_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', green_visual./255);
hold on;

[r, pval] = corr(P.vis_prefDir_diff, rtoZ(P.rsc_vis), 'Type', 'Spearman');

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_sac(P.sac_prefDir_diff>= bin_edges(b) & P.sac_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', purple_saccadey./255);
hold on;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_pur(P.pur_prefDir_diff >= bin_edges(b) & P.pur_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'go:', 'LineWidth', 2, 'Color', red_pursuity./255);

xticks(bin_centers)
xlim([5,175])
ylim([0,0.025])

prettyFig;
axis square

ax1(2) = nexttile;
P = pairs_tbl_uni(pairs_tbl_uni.n1_probe~=pairs_tbl_uni.n2_probe,:);

% even width 
%bin_edges = linspace(,7000,N_BINS+1);
bin_edges = 0:20:180;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_vis(P.vis_prefDir_diff >= bin_edges(b) & P.vis_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', green_visual./255);
hold on;

[r, pval] = corr(P.vis_prefDir_diff, rtoZ(P.rsc_vis), 'Type', 'Spearman');

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_sac(P.sac_prefDir_diff>= bin_edges(b) & P.sac_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'o-', 'LineWidth', 2, 'Color', purple_saccadey./255);
hold on;

[mean_rho,sem_rho] = deal(zeros(1,numel(bin_centers)));
for b = 1:length(bin_centers)
    rho_vals = rtoZ(P.rsc_pur(P.pur_prefDir_diff >= bin_edges(b) & P.pur_prefDir_diff < bin_edges(b+1)));
    n_pairs = numel(rho_vals);

    mean_rho(b) = mean(rho_vals, 'omitnan');
    sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
end
errorbar(bin_centers, mean_rho, sem_rho, 'go:', 'LineWidth', 2, 'Color', red_pursuity./255);

xticks(bin_centers)
xlim([5,175])
ylim([0,0.013])

prettyFig;
axis square

linkaxes(ax1,'x')

savebigPDF(f4f, '/Users/kendranoneman/Milestones/proposal/figs/rsc_tuning_crosshemi.pdf') 

%% Heatmaps of rsc, VMI quantiles (across hemi)
P = pairs_tbl(pairs_tbl.n1_probe~=pairs_tbl.n2_probe,:);

f4g = figure;
f4g.Position = [100 100 800 400]; 
tl = tiledlayout(1,3); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'compact'; 

xx = P.n1_VMI;
yy = P.n2_VMI;

% VISUAL
ax1(1) = nexttile;
cc = P.rsc_vis;

% Compute mean values per (xx, yy)
M1 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);

% Compute counts per (xx, yy)
N1 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M1);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(green_visual,0,256);
clim([0,0.0118]);
colormap(ax1(1), cmap);
cb = colorbar(ax1(1));

prettyFig;

% SACCADE
ax1(2) = nexttile;
cc = P.rsc_sac;

M2 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N2 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M2);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(purple_saccadey,0,256);
clim([0,0.0118]);
colormap(ax1(2), cmap);
cb = colorbar(ax1(2));

prettyFig;

% VISUAL
ax1(3) = nexttile;
cc = P.rsc_pur;

M3 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N3 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M3);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(red_pursuity,0,256);
clim([0,0.0118]);
colormap(ax1(3), cmap);
cb = colorbar(ax1(3));

prettyFig;

savebigPDF(f4g, '/Users/kendranoneman/Milestones/proposal/figs/rsc_heatmap_vmi_crosshemi.pdf') 

%% Heatmaps of rsc, SPI quantiles (across hemi)
P = pairs_tbl(pairs_tbl.n1_probe~=pairs_tbl.n2_probe,:);

f4h = figure;
f4h.Position = [100 100 800 400]; 
tl = tiledlayout(1,3); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'compact'; 

xx = P.n1_SPI;
yy = P.n2_SPI;

% VISUAL
ax1(1) = nexttile;
cc = P.rsc_vis;

% Compute mean values per (xx, yy)
M4 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);

% Compute counts per (xx, yy)
N4 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M4);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(green_visual,0,256);
clim([0,0.0118]);
colormap(ax1(1), cmap);
cb = colorbar(ax1(1));

prettyFig;

% SACCADE
ax1(2) = nexttile;
cc = P.rsc_sac;

M5 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N5 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M5);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(purple_saccadey,0,256);
clim([0,0.0118]);
colormap(ax1(2), cmap);
cb = colorbar(ax1(2));

prettyFig;

% VISUAL
ax1(3) = nexttile;
cc = P.rsc_pur;

M6 = accumarray([yy(:), xx(:)], cc(:), [5 5], @mean, NaN);
N6 = accumarray([yy(:), xx(:)], 1, [5 5], @sum, 0);

% Plot the heatmap
imagesc(1:5, 1:5, M6);
axis equal tight;
set(gca, 'YDir', 'normal'); % y increases upward
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
cmap = customColormap(red_pursuity,0,256);
clim([0,0.0118]);
colormap(ax1(3), cmap);
cb = colorbar(ax1(3));

prettyFig;
savebigPDF(f4h, '/Users/kendranoneman/Milestones/proposal/figs/rsc_heatmap_spi_crosshemi.pdf') 

%% Pop metrics

f4j = figure;
f4j.Position = [100 100 1500 500];
tl = tiledlayout(2,5);
tl.TileSpacing = 'loose';
tl.Padding = 'loose';

for m = 1:size(all_pop_mats,4)
    ax1(m) = nexttile;
    hold on;
    for p = 1:size(all_pop_mats,2)
        values = [all_pop_mats(:,p,1,m); all_pop_mats(:,p,2,m)];

        if p==1
            errorbar(p, mean(values,'omitnan'), std(values,'omitnan'), 'o', 'LineWidth', 2, 'Color', green_visual./255);
        elseif p==2
            errorbar(p, mean(values,'omitnan'), std(values,'omitnan'), 'o', 'LineWidth', 2, 'Color', purple_saccadey./255);
        else
            errorbar(p, mean(values,'omitnan'), std(values,'omitnan'), 'o', 'LineWidth', 2, 'Color', red_pursuity./255);
        end
        xlim([0,4])
        xticks([1,2,3])
        prettyFig;       
    end

    % fprintf('metric = %d\n', m);
    % values1 = [all_pop_mats(:,2,1,m); all_pop_mats(:,2,2,m)];
    % values2 = [all_pop_mats(:,3,1,m); all_pop_mats(:,3,2,m)];
    % 
    % nPerm = 10000;  % number of permutations
    % obs_diff = mean(values1,'omitnan') - mean(values2,'omitnan');  % observed difference in means
    % 
    % % Combine all values
    % all_vals = [values1(:); values2(:)];
    % n1 = numel(values1);
    % 
    % % Preallocate for speed
    % perm_diffs = zeros(nPerm,1);
    % 
    % for i = 1:nPerm
    %     shuffled = all_vals(randperm(numel(all_vals)));
    %     perm_diffs(i) = mean(shuffled(1:n1),'omitnan') - mean(shuffled(n1+1:end),'omitnan');
    % end
    % 
    % % Two-sided p-value
    % p = mean(abs(perm_diffs) >= abs(obs_diff),'omitnan');
    % 
    % fprintf('Observed mean difference = %.4f, p = %.4f\n', obs_diff, p);
    blah = 1;

end

for m = 1:size(all_pop_mats,4)
    ax1(m+5) = nexttile;
    hold on;
    for p = 1:size(all_pop_mats,2)
        values = all_pop_mats(:,p,3,m);

        if p==1
            errorbar(p, mean(values,'omitnan'), std(values,'omitnan'), 'o', 'LineWidth', 2, 'Color', green_visual./255);
        elseif p==2
            errorbar(p, mean(values,'omitnan'), std(values,'omitnan'), 'o', 'LineWidth', 2, 'Color', purple_saccadey./255);
        else
            errorbar(p, mean(values,'omitnan'), std(values,'omitnan'), 'o', 'LineWidth', 2, 'Color', red_pursuity./255);
        end
        xlim([0,4])
        xticks([1,2,3])
        prettyFig;       
    end

    fprintf('metric = %d\n', m);
    values1 = all_pop_mats(:,2,3,m);
    values2 = all_pop_mats(:,3,3,m);

    nPerm = 10000;  % number of permutations
    obs_diff = mean(values1,'omitnan') - mean(values2,'omitnan');  % observed difference in means
    
    % Combine all values
    all_vals = [values1(:); values2(:)];
    n1 = numel(values1);
    
    % Preallocate for speed
    perm_diffs = zeros(nPerm,1);
    
    for i = 1:nPerm
        shuffled = all_vals(randperm(numel(all_vals)));
        perm_diffs(i) = mean(shuffled(1:n1),'omitnan') - mean(shuffled(n1+1:end),'omitnan');
    end
    
    % Two-sided p-value
    p = mean(abs(perm_diffs) >= abs(obs_diff),'omitnan');
    
    fprintf('Observed mean difference = %.4f, p = %.4f\n', obs_diff, p);
end

linkaxes(ax1([1,6]),'xy')
linkaxes(ax1([2,7]),'xy')
linkaxes(ax1([3,8]),'xy')
linkaxes(ax1([4,9]),'xy')
linkaxes(ax1([5,10]),'xy')

%savebigPDF(f4j, '/Users/kendranoneman/Milestones/proposal/figs/pop_metrics.pdf') 

%% CHAPTER 5 

%% Pie charts of pure v. saccadey pursuit

TT = Tpurs(Tpurs.monkey=='scrappy' & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0),:);

np1 = 0;
for t=1:height(TT)
    [~,rh] = cart2pol(TT.eyeVel{t}(1,:),TT.eyeVel{t}(2,:));
    if TT.pursType(t)=="pure"% && max(rh(TT.pursuitOnset(t)-100:TT.pursuitOnset(t)+100))<50
        np1 = np1 + 1;
    end
end

nc1 = 0;
for t=1:height(TT)
    [~,rh] = cart2pol(TT.eyeVel{t}(1,:),TT.eyeVel{t}(2,:));
    if TT.pursType(t)~="pure"
        nc1 = nc1 + 1;
    end
end

TT = Tpurs(Tpurs.monkey=='walter' & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0),:);

perPure1 = (np1/(nc1+np1))*100;
perCS1 = 100-(np1/(nc1+np1))*100;

np2 = 0;
for t=1:height(TT)
    [~,rh] = cart2pol(TT.eyeVel{t}(1,:),TT.eyeVel{t}(2,:));
    if TT.pursType(t)=="pure"% && max(rh(TT.pursuitOnset(t)-100:TT.pursuitOnset(t)+100))<50
        np2 = np2 + 1;
    end
end

nc2 = 0;
for t=1:height(TT)
    [~,rh] = cart2pol(TT.eyeVel{t}(1,:),TT.eyeVel{t}(2,:));
    if TT.pursType(t)~="pure"
        nc2 = nc2 + 1;
    end
end

perPure2 = (np2/(nc2+np2))*100;
perCS2 = 100-(np2/(nc2+np2))*100;

%%
f5aa = figure;
f5aa.Position = [100 100 500 300];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

ax1(1) = nexttile;
values = [perPure1, perCS1];
pie(values);
colormap([0 0 0; 0.3490 0.3490 0.3490]); % optional custom colors
prettyFig;

ax1(2) = nexttile;
values = [perPure2, perCS2];
pie(values);
colormap([0 0 0; 0.3490 0.3490 0.3490]); % optional custom colors
prettyFig;

savebigPDF(f5aa, '/Users/kendranoneman/Milestones/proposal/figs/cs_pies.pdf')

%% Behavior traces
TT = Tpurs(Tpurs.sess_name=='kendra_scrappy_0140a' & (isnan(Tpurs.msOffset) | Tpurs.msOffset<-100),:);

x = -200:1:200;

f5a = figure;
f5a.Position = [100 100 1100 200]; 
tl = tiledlayout(1,2); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'compact';

ax1(1) = nexttile;
np = 0;
for t=1:height(TT)
    [~,rh] = cart2pol(TT.eyeVel{t}(1,:),TT.eyeVel{t}(2,:));

    if TT.pursType(t)=="pure" && max(rh(TT.pursuitOnset(t)-200:TT.pursuitOnset(t)+200))<50
        plot(x,rh(TT.pursuitOnset(t)-200:TT.pursuitOnset(t)+200),'k-')
        np = np + 1;
    end
    hold on;
end
ylim([0,200])
prettyFig;

ax1(2) = nexttile;
nc = 0;
for t=1:height(TT)
    [~,rh] = cart2pol(TT.eyeVel{t}(1,:),TT.eyeVel{t}(2,:));

    if TT.pursType(t)~="pure"
        plot(x,rh(TT.pursuitOnset(t)-200:TT.pursuitOnset(t)+200),'k-')
        nc = nc + 1;
    end
    hold on;
end
prettyFig;
ylim([0,200])

linkaxes(ax1,'xy')

savebigPDF(f5a, '/Users/kendranoneman/Milestones/proposal/figs/catchup_behav.pdf')

%% PSTH colored by VMI and SPI quantiles (for pure v. cs pursuit)
nBins = 5;
sigma = 10;
visWindow = [-100,300];
sacWindow = [-200,200];
purWindow = [-200,200];

CC.VMI_quantile = discretize(CC.VMI,  quantile(CC.VMI, linspace(0, 1, nBins+1)));
CC.VMI_quantile(CC.VMI == max(CC.VMI)) = nBins;

[vis_mn_vmi,vis_sem_vmi,sac_mn_vmi,sac_sem_vmi,pur_mn_vmi,pur_sem_vmi] = deal(cell(18,nBins));
for q = 1:nBins
    fprintf('Bin %d of %d...\n', q, nBins);
    these_units = CC(CC.VMI_quantile==q,:);

    these_sess = unique(these_units.sess_name);
    for s = 1:numel(these_sess)
        this_pursTbl1 = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
        this_pursTbl2 = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType~='pure',:);
        
        % PURS PSTH
        time_window = purWindow;

        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl1.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl1.spiketimes_2, 'uni', 0);

        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_pursTbl1.pursuitOnset), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_pursTbl1.pursuitOnset), 'uni', 0);

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

        % SACC PSTH
        time_window = purWindow;

        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl2.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl2.spiketimes_2, 'uni', 0);

        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_pursTbl2.pursuitOnset), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_pursTbl2.pursuitOnset), 'uni', 0);

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

        sac_mn_vmi{s,q} = pur_mn;
        sac_sem_vmi{s,q} = pur_sem;

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
        this_pursTbl1 = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);
        this_pursTbl2 = Tpurs(Tpurs.sess_name==these_sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType~='pure',:);

        % PURS PSTH
        time_window = purWindow;

        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl1.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl1.spiketimes_2, 'uni', 0);

        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_pursTbl1.pursuitOnset), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_pursTbl1.pursuitOnset), 'uni', 0);

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

        % SAC PSTH
        time_window = purWindow;

        these_spikes1 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==1)+1), this_pursTbl2.spiketimes_1, 'uni', 0);
        these_spikes2 = cellfun(@(q) q(these_units.cluster_id(these_units.sess_name==these_sess(s) & these_units.probe_index==2)+1), this_pursTbl2.spiketimes_2, 'uni', 0);

        spike_times1 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes1, num2cell(this_pursTbl2.pursuitOnset), 'uni', 0);
        spike_times2 = cellfun(@(w,v) cellfun(@(q) q-v, w, 'uni', 0), these_spikes2, num2cell(this_pursTbl2.pursuitOnset), 'uni', 0);

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

        sac_mn_spi{s,q} = pur_mn;
        sac_sem_spi{s,q} = pur_sem;
    end
end

%%
vmi_cmap = customColormap(purple_saccadey, green_visual, nBins);

[psth_mn,psth_sem] = deal(zeros(nBins,800));
for q = 1:nBins
    psth_mn_all = [vertcat(pur_mn_vmi{:,q}) vertcat(sac_mn_vmi{:,q})];
    %psth_mn_norm = psth_mn_all;
    psth_mn_norm = zscore(psth_mn_all,[],2);
    %psth_mn_norm = psth_mn_all ./ max(psth_mn_all, [], 2);

    psth_mn(q,:) = mean(psth_mn_norm,1);
    
    psth_sem_all = [vertcat(pur_sem_vmi{:,q}) vertcat(sac_sem_vmi{:,q})];
    %psth_sem_norm = psth_sem_all;
    psth_sem_norm = psth_sem_all ./ (2*std(psth_mn_all, 0, 2));  % divide by the same std as z-scoring
    %psth_sem_norm = psth_sem_all ./ max(psth_mn_all, [], 2);

    psth_sem(q,:) = mean(psth_sem_norm,1);
end

f5b = figure;
f5b.Position = [100 100 1100 300];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

ax1(1) = nexttile;
x = purWindow(1):(purWindow(2)-1);

hold on;
for q = 1:nBins
    y = psth_mn(q,1:400);
    yu = psth_mn(q,1:400)+psth_sem(q,1:400); yl = psth_mn(q,1:400)-psth_sem(q,1:400);

    fill([x fliplr(x)], [yu fliplr(yl)], vmi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',vmi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

ax1(2) = nexttile;
x = purWindow(1):(purWindow(2)-1);

hold on;
for q = 1:nBins
    y = psth_mn(q,401:800);
    yu = psth_mn(q,401:800)+psth_sem(q,401:800); yl = psth_mn(q,401:800)-psth_sem(q,401:800);

    fill([x fliplr(x)], [yu fliplr(yl)], vmi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',vmi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

linkaxes(ax1,'y')

savebigPDF(f5b, '/Users/kendranoneman/Milestones/proposal/figs/psth_vmi_catchup.pdf')

%%
spi_cmap = customColormap(red_pursuity,purple_saccadey, nBins);

[psth_mn,psth_sem] = deal(zeros(nBins,800));
for q = 1:nBins
    psth_mn_all = [vertcat(pur_mn_spi{:,q}) vertcat(sac_mn_spi{:,q})];
    psth_mn_norm = zscore(psth_mn_all,[],2);

    psth_mn(q,:) = mean(psth_mn_norm,1);
    
    psth_sem_all = [vertcat(pur_sem_spi{:,q}) vertcat(sac_sem_spi{:,q})];
    psth_sem_norm = psth_sem_all ./ (2*std(psth_mn_all, 0, 2)); % divide by the same std as z-scoring

    psth_sem(q,:) = mean(psth_sem_norm,1);
end

f5c = figure;
f5c.Position = [100 100 1100 300];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

ax1(1) = nexttile;
x = purWindow(1):(purWindow(2)-1);

hold on;
for q = 1:nBins
    y = psth_mn(q,1:400);
    yu = psth_mn(q,1:400)+psth_sem(q,1:400); yl = psth_mn(q,1:400)-psth_sem(q,1:400);

    fill([x fliplr(x)], [yu fliplr(yl)], spi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',spi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

ax1(2) = nexttile;
x = purWindow(1):(purWindow(2)-1);

hold on;
for q = 1:nBins
    y = psth_mn(q,401:800);
    yu = psth_mn(q,401:800)+psth_sem(q,401:800); yl = psth_mn(q,401:800)-psth_sem(q,401:800);

    fill([x fliplr(x)], [yu fliplr(yl)], spi_cmap(q,:), 'linestyle', 'none', 'FaceAlpha', 0.1);
    plot(x,y,'-','Color',spi_cmap(q,:),'LineWidth',3);
    
end
prettyFig;

linkaxes(ax1,'y')

savebigPDF(f5c, '/Users/kendranoneman/Milestones/proposal/figs/psth_spi_catchup.pdf')


%% OTHER PLOTS NOT USING RN

%% Distributions of preferred directions (per monkey & probe)
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


%% Rsc and eye movement direction
P = pairs_tbl(pairs_tbl.n1_probe==pairs_tbl.n2_probe,:);

f3i = figure;
f3i.Position = [100 100 1200 400]; 
tl = tiledlayout(2,3); 
tl.TileSpacing = 'compact'; 
tl.Padding = 'compact'; 

% VISUAL
ax1(1) = nexttile(1);

xx = P.fr_perDir_vis;
yy = discretize(P.vis_prefDir_diff, 0:20:180);

[M1,M1_sem] = deal(nan(numel(unique(yy)),8));
for y = 1:numel(unique(yy))
    M1(y,:) = mean(xx(yy==y,:),'omitnan');
    M1_sem(y,:) = std(xx(yy==y,:), 0, 1, 'omitnan') ./ sqrt(sum(~isnan(xx(yy==y,:)), 1));
end


% Plot the heatmap
imagesc(-180:45:135, 10:20:170, M1);
axis tight;
%set(gca, 'YDir', 'normal'); % y increases upward
%set(gca, 'XTickLabel', []);
%set(gca, 'YTickLabel', []);
cmap = customColormap(green_visual,0,256);
%clim([0.0036,0.0248]);
colormap(ax1(1), cmap);
cb = colorbar(ax1(1));

xticks([-180,-90,0,90])
yticks([10,50,90,130,170])

prettyFig;

ax2(1) = nexttile(4);
errorbar(-180:45:135, mean(M1), mean(M1_sem),'ko-')
prettyFig;


% SACCADE
ax1(2) = nexttile(2);

xx = P.fr_perDir_sac;
yy = discretize(P.sac_prefDir_diff, 0:20:180);

[M2,M2_sem] = deal(nan(numel(unique(yy)),8));
for y = 1:numel(unique(yy))
    M2(y,:) = mean(xx(yy==y,:),'omitnan');
    M2_sem(y,:) = std(xx(yy==y,:), 0, 1, 'omitnan') ./ sqrt(sum(~isnan(xx(yy==y,:)), 1));
end

% Plot the heatmap
imagesc(-180:45:135, 10:20:170, M2);
axis tight;
%set(gca, 'YDir', 'normal'); % y increases upward
%set(gca, 'XTickLabel', []);
%set(gca, 'YTickLabel', []);
cmap = customColormap(purple_saccadey,0,256);
%clim([0.0036,0.0248]);
colormap(ax1(2), cmap);
cb = colorbar(ax1(2));

xticks([-180,-90,0,90])
yticks([10,50,90,130,170])

prettyFig;

ax2(1) = nexttile(5);
errorbar(-180:45:135, mean(M2), mean(M2_sem),'ko-')
prettyFig;

% PURSUIT
ax1(3) = nexttile(3);

xx = P.fr_perDir_pur;
yy = discretize(P.pur_prefDir_diff, 0:20:180);

[M3,M3_sem] = deal(nan(numel(unique(yy)),8));
for y = 1:numel(unique(yy))
    M3(y,:) = mean(xx(yy==y,:),'omitnan');
    M3_sem(y,:) = std(xx(yy==y,:), 0, 1, 'omitnan') ./ sqrt(sum(~isnan(xx(yy==y,:)), 1));
end

% Plot the heatmap
imagesc(-180:45:135, 10:20:170, M3);
axis tight;
%set(gca, 'YDir', 'normal'); % y increases upward
%set(gca, 'XTickLabel', []);
%set(gca, 'YTickLabel', []);
cmap = customColormap(red_pursuity,0,256);
%clim([0.0036,0.0248]);
colormap(ax1(3), cmap);
cb = colorbar(ax1(3));

xticks([-180,-90,0,90])
yticks([10,50,90,130,170])

prettyFig;

ax2(1) = nexttile(6);
errorbar(-180:45:135, mean(M3), mean(M3_sem),'ko-')
prettyFig;





