%% GRC FIGURES
data_path = '/Users/kendranoneman/Data/dualhemi_unleashed';

%% individual session

sess = 'kendra_scrappy_0142a_g0';
load(fullfile(data_path,[sess '_unleashed.mat']), 'S')

EXAMPLE_CLUSTERS = {[0,366],[0,1],[0,295],[0,325]}; % [imec, cluster_id] 

%sv = 44, 73
%sm = 36, 65, 70, 119, 257, 264, 283, 295*, 325*, 327, 335, 352
%pv = 1*, 85, 94, 224, 261, 270, 279, 280
%pm = 59, 234*, 239*, 271, 294*, 366*

%S.sessionName = sess;

%% adding some additional thresholds
FR_thresh     =  1;    % Hz
seldir_range =  [0 0.99]; 
msOff_thresh  =  -50; % ms
sacLat_range =  [100 300];  % ms

%---- Behavior criteria ----%
% mdir criteria 
Tmdir = S.mdir1.tbl(S.mdir1.tbl.result=='CORRECT',:);
Tmdir.saccadeLatency = Tmdir.SACCADE-cellfun(@(q) q(1), Tmdir.FIX_OFF);
Tmdir = Tmdir(Tmdir.saccadeLatency >= sacLat_range(1) & Tmdir.saccadeLatency < sacLat_range(2),:);

[th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tmdir.eyePos, 'uni', 0);
Tmdir = Tmdir(cellfun(@(r,t,s,a) (max(r(s-10:s+75)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s-10:s+75))) - a)) < 45), rh, th, num2cell(Tmdir.SACCADE), num2cell(Tmdir.angle), 'uni', 1),:);

% purs criteria
Tpurs = S.purs1.tbl(S.purs1.tbl.result=='CORRECT' & S.purs1.tbl.pursType=='pure',:);
Tpurs = Tpurs(isnan(Tpurs.msOffset) | Tpurs.msOffset<msOff_thresh,:);

[th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tpurs.eyePos, 'uni', 0);
Tpurs = Tpurs(cellfun(@(r,t,s,a) (max(r(s+200:s+800)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s+200:s+800))) - a)) < 160), rh, th, num2cell(Tpurs.pursuitOnset), num2cell(Tpurs.angle), 'uni', 1),:);

%---- Cluster criteria ----%
Cleft = S.kilosort(1).clusters; Crght = S.kilosort(2).clusters;

% firing rate
Cleft = Cleft(Cleft.mdir1_Hz > FR_thresh & Cleft.purs1_Hz > FR_thresh,:);
Crght = Crght(Crght.mdir1_Hz > FR_thresh & Crght.purs1_Hz > FR_thresh,:);

% selectivity
Cleft = Cleft((Cleft.vis_sel_dir>seldir_range(1) & Cleft.vis_sel_dir<seldir_range(2)) & (Cleft.sac_sel_dir>seldir_range(1) & Cleft.sac_sel_dir<seldir_range(2)) & (Cleft.pur_sel_dir>seldir_range(1) & Cleft.pur_sel_dir<seldir_range(2)),:);
Crght = Crght((Crght.vis_sel_dir>seldir_range(1) & Crght.vis_sel_dir<seldir_range(2)) & (Crght.sac_sel_dir>seldir_range(1) & Crght.sac_sel_dir<seldir_range(2)) & (Crght.pur_sel_dir>seldir_range(1) & Crght.pur_sel_dir<seldir_range(2)),:);

%% FIGURE 2: Plotting waveform templates 

f2 = figure;
f2.Position = [100 100 300 1200];

for c = 1:height(Cleft)
    cinds = S.kilosort(1).spike_clusters == Cleft.cluster_id(c);
    temps = unique(S.kilosort(1).spike_templates(cinds));

    mean_wf = mean(S.kilosort(1).templates(temps+1,:,:),3);
    mean_wf = mean_wf(1:end-10);
    spike_pos = mean(S.kilosort(1).spike_positions(cinds,:),1);

    x_center = (spike_pos(1) - 11 - (103/2))./3;

    x = ((0:(numel(mean_wf)-1))/10) + x_center;
    x = (x + randi([1 20])./2);
    y = (mean_wf - mean_wf(1)).*1.5;
    y = y + (spike_pos(2)./1000);
    %y = (mean_wf + spike_pos(2)/1000);


    plot(x,y,'k-','linewidth',0.6)
    hold on;

    blah = 1;

end
ylim([-0.3420 (max(S.kilosort(1).spike_positions(:,2))./1000)+0.2]);
prettyFig;

savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/waveforms_probe.pdf')

%% FIGURE 3: Plot polar eye traces for tasks
PREINT = 10; POSTINT = 75;

f3a = figure;
T = Tmdir(Tmdir.distance==10,:);

% Loop through each trial in the table
for t = 1:height(T)
    this_angle    =  T.angle(t);   % angle in degrees
    this_saccTime =  T.SACCADE(t); % time (ms) of saccade

    % Pull out eye traces for trial and convert to polar coordinates
    this_eyePos = T.eyePos{t}; % HE, VE eye traces
    [theta,rho] = cart2pol(T.eyePos{t}(1,:),T.eyePos{t}(2,:));

    % Adding a check here to not include trials where his first saccade was in the wrong direction or way beyond the bounds of the screen
    polarplot(theta(this_saccTime-PREINT:this_saccTime+POSTINT),rho(this_saccTime-PREINT:this_saccTime+POSTINT), ...
            'color',[0,0,0],'linewidth',1)
    hold on;
end
rlim([0 25]);
thetaticks(0:45:315);
prettyFig;
savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/sacc_10deg.pdf')

%--------------------------------------
f3b = figure;
T = Tmdir(Tmdir.distance==20,:);

% Loop through each trial in the table
for t = 1:height(T)
    this_angle    =  T.angle(t);   % angle in degrees
    this_saccTime =  T.SACCADE(t); % time (ms) of saccade

    % Pull out eye traces for trial and convert to polar coordinates
    this_eyePos = T.eyePos{t}; % HE, VE eye traces
    [theta,rho] = cart2pol(T.eyePos{t}(1,:),T.eyePos{t}(2,:));

    % Adding a check here to not include trials where his first saccade was in the wrong direction or way beyond the bounds of the screen
    
    polarplot(theta(this_saccTime-PREINT:this_saccTime+POSTINT),rho(this_saccTime-PREINT:this_saccTime+POSTINT), ...
        'color',[0,0,0],'linewidth',1)
    hold on;
   
end
rlim([0 25]);
thetaticks(0:45:315);
prettyFig;
savebigPDF(2, '/Users/kendranoneman/Posters/GRC-2025/sacc_20deg.pdf')

%--------------------------------------
%--------------------------------------
PREINT = 0; POSTINT = 1000;
f3c = figure;
T = Tpurs(Tpurs.pursuitSpeed==15,:);

% Loop through each trial in the table
for t = 1:height(T)
    this_angle    =  T.angle(t);   % angle in degrees
    this_pursTime =  T.pursuitOnset(t); % time (ms) of saccade

    % Pull out eye traces for trial and convert to polar coordinates
    this_eyePos = T.eyePos{t}; % HE, VE eye traces
    [theta,rho] = cart2pol(T.eyePos{t}(1,:),T.eyePos{t}(2,:));

    % Adding a check here to not include trials where his first saccade was in the wrong direction or way beyond the bounds of the screen
    polarplot(theta(this_pursTime-PREINT:this_pursTime+POSTINT),rho(this_pursTime-PREINT:this_pursTime+POSTINT), ...
            'color',[0,0,0],'linewidth',1)
    hold on;
end
rlim([0 25]);
thetaticks(0:45:315);
prettyFig;
savebigPDF(3, '/Users/kendranoneman/Posters/GRC-2025/purs_15deg.pdf')

%--------------------------------------
f3d = figure;
T = Tpurs(Tpurs.pursuitSpeed==20,:);

% Loop through each trial in the table
for t = 1:height(T)
    this_angle    =  T.angle(t);   % angle in degrees
    this_pursTime =  T.pursuitOnset(t); % time (ms) of saccade

    % Pull out eye traces for trial and convert to polar coordinates
    this_eyePos = T.eyePos{t}; % HE, VE eye traces
    [theta,rho] = cart2pol(T.eyePos{t}(1,:),T.eyePos{t}(2,:));

    % Adding a check here to not include trials where his first saccade was in the wrong direction or way beyond the bounds of the screen
    polarplot(theta(this_pursTime-PREINT:this_pursTime+POSTINT),rho(this_pursTime-PREINT:this_pursTime+POSTINT), ...
            'color',[0,0,0],'linewidth',1)
    hold on;
end
rlim([0 25]);
thetaticks(0:45:315);
prettyFig;
savebigPDF(4, '/Users/kendranoneman/Posters/GRC-2025/purs_20deg.pdf')

%% FIGURE 3: Plotting eye velocity  to point out epochs

f4a = figure;
f4a.Position = [100 100 900 300];

%tl = tiledlayout(2,2);
%tl.TileSpacing = 'compact';
%tl.Padding = 'tight';

% mdir, aligned to stim onset
%nexttile
x = -300:1:500;
FR_WIN = [50 150];
fill([FR_WIN fliplr(FR_WIN)], [[1000 1000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.15);
hold on;
xline(0,'k--')

for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyeVel{t}(1,:),Tmdir.eyeVel{t}(2,:));
    plot(x,rh(Tmdir.TARG_ON{t}(1)-300:Tmdir.TARG_ON{t}(1)+500),'k-')
end
%xlabel('time aligned to stimulus onset (ms)')
%ylabel('radial eye velocity (deg/s)')
ylim([0 1000])
prettyFig;
savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/sacc_eyes_stim.pdf')

ia_mdirRasters(S, 'IMEC', 0, 'ALIGN', 'stim', 'CLUSTER', 325, 'TICK_LENGTH', 2, 'Y_LIMITS', [0 200])
savebigPDF(2, '/Users/kendranoneman/Posters/GRC-2025/sacc_spks_stim.pdf')


f4b = figure;
f4b.Position = [100 100 900 300];

x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[1000 1000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.15);
hold on;
xline(0,'k--')
for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyeVel{t}(1,:),Tmdir.eyeVel{t}(2,:));
    plot(x,rh(Tmdir.SACCADE(t)-300:Tmdir.SACCADE(t)+500),'k-')
end
%xlabel('time aligned to saccade onset (ms)')
%ylabel('radial eye velocity (deg/s)')
ylim([0 1000])
%title('motor epoch (MGS)')
prettyFig;
savebigPDF(3, '/Users/kendranoneman/Posters/GRC-2025/sacc_eyes_sacc.pdf')

ia_mdirRasters(S, 'IMEC', 0, 'ALIGN', 'sacc', 'CLUSTER', 325, 'TICK_LENGTH', 2, 'Y_LIMITS', [0 200])
savebigPDF(4, '/Users/kendranoneman/Posters/GRC-2025/sacc_spks_sacc.pdf')

f4c = figure;
f4c.Position = [100 100 900 300];

x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[100 100] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tpurs)
    [~,rh] = cart2pol(Tpurs.eyeVel{t}(1,:),Tpurs.eyeVel{t}(2,:));
    [pursuit_onset, rxnTime, msOffset, csOnset, csVelocity, csPeak, csOffset, csAngle, csType] = detect_pursuitOnset(Tpurs.eyePos{t}, Tpurs.eyeVel{t}, Tpurs.PURSUIT_TARG_ON(t), S.purs1.params.crossingTime, Tpurs.pursuitSpeed(t), Tpurs.angle(t), 'CS_PREINT', 50, 'CS_POSTINT', 100);
    if isequal(csType,'pure')
        plot(x,rh(Tpurs.pursuitOnset(t)-300:Tpurs.pursuitOnset(t)+500),'k-')
    end
    Tpurs.pursuitOnset(t) = pursuit_onset;
    Tpurs.pursuitLatency(t) = rxnTime;
    Tpurs.pursType(t) = categorical(string(csType));

end
%xlabel('time aligned to pursuit onset (ms)')
%ylabel('radial eye velocity (deg/s)')
ylim([0 50])
%title('motor epoch (pursuit)')
prettyFig;
savebigPDF(5, '/Users/kendranoneman/Posters/GRC-2025/purs_eyes.pdf')


ia_pursRasters(S, 'IMEC', 0, 'ALIGN', 'purs', 'PURE_ONLY', true, 'CLUSTER', 325, 'TICK_LENGTH', 2, 'Y_LIMITS', [0 200])
savebigPDF(6, '/Users/kendranoneman/Posters/GRC-2025/purs_spks.pdf')

%% FIGURE 4: Example units 

clust_id = 366;

ia_mdirRasters(S, 'IMEC', 0, 'ALIGN', 'stim', 'CLUSTER', clust_id, 'TICK_LENGTH', 2, 'X_LIMITS', [-200 200]);%, 'Y_LIMITS', [0 110]);
%savebigPDF(1, sprintf('/Users/kendranoneman/Posters/GRC-2025/clust%.4d_stim.pdf',clust_id))

ia_mdirRasters(S, 'IMEC', 0, 'ALIGN', 'sacc', 'CLUSTER', clust_id, 'TICK_LENGTH', 2, 'X_LIMITS', [-200 200]);%, 'Y_LIMITS', [0 110]);
%savebigPDF(2, sprintf('/Users/kendranoneman/Posters/GRC-2025/clust%.4d_sacc.pdf',clust_id))

ia_pursRasters(S, 'IMEC', 0, 'ALIGN', 'purs', 'PURE_ONLY', true, 'CLUSTER', clust_id, 'TICK_LENGTH', 2, 'X_LIMITS', [-200 200]);%, 'Y_LIMITS', [0 110]);
%savebigPDF(3, sprintf('/Users/kendranoneman/Posters/GRC-2025/clust%.4d_purs.pdf',clust_id))


%% compiling all of the sessions together

files = dir(fullfile(data_path, '*.mat'));
files = files(~strcmp({files.name}, 'COMBINED.mat'));

if isfile(fullfile(data_path, 'COMBINED.mat'))
    load(fullfile(data_path, 'COMBINED.mat'), 'CLUSTERS', 'MDIR', 'PURS');
    fprintf('Loaded existing COMBINED.mat\n');
else
    fprintf('Creating COMBINED.mat\n');
    CLUSTERS = table();
    MDIR = struct(); PURS = struct();
    
    for i = 1:length(files)
        file_name = files(i).name;
        file_path = fullfile(data_path, file_name);

        fprintf('Processing file: %s\n', file_name); 
    
        load(file_path, 'S');  % Load only the variable S
    
        ks = S.kilosort;
    
        % --- clusters --- %
        for k = 1:numel(ks)
            if isfield(ks(k), 'clusters') && istable(ks(k).clusters)

                new_tbl = ks(k).clusters;
               
                if i>1 || k>1
                    missing_in_clusters = setdiff(new_tbl.Properties.VariableNames, CLUSTERS.Properties.VariableNames);
                    for v = missing_in_clusters
                        CLUSTERS.(v{1}) = NaN(height(CLUSTERS),1);
                    end
                    
                    % Add any vars in CLUSTERS that new_tbl is missing
                    missing_in_new = setdiff(CLUSTERS.Properties.VariableNames, new_tbl.Properties.VariableNames);
                    for v = missing_in_new
                        new_tbl.(v{1}) = NaN(height(new_tbl),1);
                    end
                    
                    % Reorder new_tbl to match CLUSTERS
                    new_tbl = new_tbl(:, CLUSTERS.Properties.VariableNames);
                end
                
                % Append
                CLUSTERS = [CLUSTERS; new_tbl];
            end
        end
    
        % --- mdir --- %
        concat_tbl = table();
        s_fields = fieldnames(S);
        for f = 1:numel(s_fields)
            fname = s_fields{f};
            if contains(fname, {'dirmem', 'mdir'})
                this_struct = S.(fname);
                if isstruct(this_struct) && isfield(this_struct, 'tbl') && istable(this_struct.tbl)
                    concat_tbl = [concat_tbl; this_struct.tbl];
                end
            end
        end
        [~, base_name, ~] = fileparts(file_name);
        base_name = erase(base_name, '_unleashed');
        if ~isempty(concat_tbl)
            MDIR.(base_name) = concat_tbl;
        end
    
        % --- purs --- %
        concat_tbl = table();
        s_fields = fieldnames(S);
        for f = 1:numel(s_fields)
            fname = s_fields{f};
            if contains(fname, {'purs', 'pursuit'})
                this_struct = S.(fname);
                if isstruct(this_struct) && isfield(this_struct, 'tbl') && istable(this_struct.tbl)
                    concat_tbl = [concat_tbl; this_struct.tbl];
                end
            end
        end
        [~, base_name, ~] = fileparts(file_name);
        base_name = erase(base_name, '_unleashed');
        if ~isempty(concat_tbl)
            PURS.(base_name) = concat_tbl;
        end
    
        fprintf('---------------\n');
    end
    
    % After the loop ends
    save(fullfile(data_path, 'COMBINED.mat'), 'CLUSTERS', 'MDIR', 'PURS', '-v7.3');
end

%% THRESHOLDS
FR_thresh     =  0.5;    % Hz
seldir_range =  [0 0.99]; 
msOff_thresh  =  -50; % ms
sacLat_range =  [100 300];  % ms

sessions = fieldnames(PURS);
CLUSTERS = CLUSTERS(ismember(CLUSTERS.sessionName,sessions),:);

TBL = [];
for f = 1:numel(sessions)
    % --- mdir --- %
    Tmdir = MDIR.(sessions{f});
    Tmdir = Tmdir(Tmdir.result=='CORRECT',:);

    Tmdir.saccadeLatency = Tmdir.SACCADE-cellfun(@(q) q(1), Tmdir.FIX_OFF);
    Tmdir = Tmdir(Tmdir.saccadeLatency >= sacLat_range(1) & Tmdir.saccadeLatency < sacLat_range(2),:);
    
    [th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tmdir.eyePos, 'uni', 0);
    Tmdir = Tmdir(cellfun(@(r,t,s,a) (max(r(s-10:s+75)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s-10:s+75))) - a)) < 45), rh, th, num2cell(Tmdir.SACCADE), num2cell(Tmdir.angle), 'uni', 1),:);

    % --- purs --- %
    Tpurs = PURS.(sessions{f});
    Tpurs = Tpurs(Tpurs.result=='CORRECT' & Tpurs.pursType=='pure',:);
    Tpurs = Tpurs(isnan(Tpurs.msOffset) | Tpurs.msOffset<msOff_thresh,:);
    
    [th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tpurs.eyePos, 'uni', 0);
    Tpurs = Tpurs(cellfun(@(r,t,s,a) (max(r(s+200:s+800)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s+200:s+800))) - a)) < 160), rh, th, num2cell(Tpurs.pursuitOnset), num2cell(Tpurs.angle), 'uni', 1),:);

    new_row = table(Tmdir.sessionName(1), {Tmdir}, {Tpurs}, 'VariableNames', {'session','mdir','purs'});
    TBL = [TBL; new_row];
end

CLUSTERS.mdir1_Hz(isnan(CLUSTERS.mdir1_Hz)) = CLUSTERS.dirmem_withhelp_varDelays_0001_Hz(isnan(CLUSTERS.mdir1_Hz));
CLUSTERS.purs1_Hz(isnan(CLUSTERS.purs1_Hz)) = CLUSTERS.pursuit_task_0001_Hz(isnan(CLUSTERS.purs1_Hz));

%% Re-calculating SPI and other metrics 
dirs = 0:45:315;

% re-calculating SPI and VMI using only "good" trials 
CLUSTS_FILT = [];
for f = 1:numel(sessions)
    this_sess = sessions{f};

    Tmdir = TBL.mdir{TBL.session==this_sess};
    Tpurs = TBL.purs{TBL.session==this_sess};
    Tclus_all = CLUSTERS(CLUSTERS.sessionName==this_sess,:);

    for i = 1:numel(unique(Tclus_all.imec))
        Tclus = Tclus_all(Tclus_all.imec==i-1,:);

        [VMIs_all,SPIs_all,VMIs_all_dp,SPIs_all_dp,VMIs_all_bc,SPIs_all_bc,VMIs_best_highEpo,SPIs_best_highEpo,VMIs_best_eachEpo,SPIs_best_eachEpo] = deal(nan(height(Tclus),1));
        [VMIs_eachDir, SPIs_eachDir] = deal(cell(height(Tclus),numel(dirs)));
        for c = 1: height(Tclus)
            % all conditions
            mdir_spks = cellfun(@(q) q{c}, Tmdir.(sprintf('spiketimes_imec%d',Tclus.imec(c))), 'uni', 0);
            purs_spks = cellfun(@(q) q{c}, Tpurs.(sprintf('spiketimes_imec%d',Tclus.imec(c))), 'uni', 0);
    
            stim_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=50 & (q-t(1))<150)/100)*1000, mdir_spks, Tmdir.TARG_ON, 'uni', 1),'omitnan');
            sacc_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, mdir_spks, num2cell(Tmdir.SACCADE), 'uni', 1),'omitnan');
            purs_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, purs_spks, num2cell(Tpurs.pursuitOnset), 'uni', 1),'omitnan');
    
            VMIs_all(c) = (stim_hz - sacc_hz) / (stim_hz + sacc_hz);
            SPIs_all(c) = (sacc_hz - purs_hz) / (sacc_hz + purs_hz);

            stim_var = var(cellfun(@(q,t) (sum((q-t(1))>=50 & (q-t(1))<150)/100)*1000, mdir_spks, Tmdir.TARG_ON, 'uni', 1),'omitnan');
            sacc_var = var(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, mdir_spks, num2cell(Tmdir.SACCADE), 'uni', 1),'omitnan');
            purs_var = var(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, purs_spks, num2cell(Tpurs.pursuitOnset), 'uni', 1),'omitnan');

            VMIs_all_dp(c) = (stim_hz - sacc_hz) / sqrt(stim_var * sacc_var);
            SPIs_all_dp(c) = (sacc_hz - purs_hz) / sqrt(sacc_var * purs_var);

            stim_base_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-400 & (q-t(1))<-200)/100)*1000, mdir_spks, Tmdir.TARG_ON, 'uni', 1),'omitnan');
            sacc_base_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-300 & (q-t(1))<-200)/100)*1000, mdir_spks, num2cell(Tmdir.SACCADE), 'uni', 1),'omitnan');
            purs_base_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-300 & (q-t(1))<-200)/100)*1000, purs_spks, num2cell(Tpurs.pursuitOnset), 'uni', 1),'omitnan');

            stim_bc_hz = max(0,stim_hz - stim_base_hz);
            sacc_bc_hz = max(0,sacc_hz - sacc_base_hz);
            purs_bc_hz = max(0,purs_hz - purs_base_hz);

            VMIs_all_bc(c) = (stim_bc_hz - sacc_bc_hz) / (stim_bc_hz + sacc_bc_hz);
            SPIs_all_bc(c) = (sacc_bc_hz - purs_bc_hz) / (sacc_bc_hz + purs_bc_hz);

            % % best direction for highest FR epoch
            % [~,max_epoch] = max([stim_hz sacc_hz purs_hz]);
            % if max_epoch==1
            %     pref_dir = Tclus.vis_pref_dir(c);
            % elseif max_epoch==2
            %     pref_dir = Tclus.sac_pref_dir(c);
            % else
            %     pref_dir = Tclus.pur_pref_dir(c);
            % end
            % [~, dir_idx] = min(abs(mod(abs(pref_dir - dirs + 180), 360) - 180));
            % best_dir = dirs(dir_idx);
            % 
            % mdir_spks = cellfun(@(q) q{c}, Tmdir.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tmdir.angle==best_dir), 'uni', 0);
            % purs_spks = cellfun(@(q) q{c}, Tpurs.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tpurs.angle==best_dir), 'uni', 0);
            % 
            % stim_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=50 & (q-t(1))<150)/100)*1000, mdir_spks, Tmdir.TARG_ON(Tmdir.angle==best_dir), 'uni', 1),'omitnan');
            % sacc_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, mdir_spks, num2cell(Tmdir.SACCADE(Tmdir.angle==best_dir)), 'uni', 1),'omitnan');
            % purs_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, purs_spks, num2cell(Tpurs.pursuitOnset(Tpurs.angle==best_dir)), 'uni', 1),'omitnan');
            % 
            % VMIs_best_highEpo(c) = (stim_hz - sacc_hz) / (stim_hz + sacc_hz);
            % SPIs_best_highEpo(c) = (sacc_hz - purs_hz) / (sacc_hz + purs_hz);
            % 
            % % best direction for each epoch
            % [~, dir_idx] = min(abs(mod(abs(Tclus.vis_pref_dir(c) - dirs + 180), 360) - 180));
            % vis_best_dir = dirs(dir_idx);
            % 
            % [~, dir_idx] = min(abs(mod(abs(Tclus.sac_pref_dir(c) - dirs + 180), 360) - 180));
            % sac_best_dir = dirs(dir_idx);
            % 
            % [~, dir_idx] = min(abs(mod(abs(Tclus.pur_pref_dir(c) - dirs + 180), 360) - 180));
            % pur_best_dir = dirs(dir_idx);
            % 
            % vis_mdir_spks = cellfun(@(q) q{c}, Tmdir.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tmdir.angle==vis_best_dir), 'uni', 0);
            % sac_mdir_spks = cellfun(@(q) q{c}, Tmdir.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tmdir.angle==sac_best_dir), 'uni', 0);
            % purs_spks = cellfun(@(q) q{c}, Tpurs.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tpurs.angle==pur_best_dir), 'uni', 0);
            % 
            % stim_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=50 & (q-t(1))<150)/100)*1000, vis_mdir_spks, Tmdir.TARG_ON(Tmdir.angle==vis_best_dir), 'uni', 1),'omitnan');
            % sacc_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, sac_mdir_spks, num2cell(Tmdir.SACCADE(Tmdir.angle==sac_best_dir)), 'uni', 1),'omitnan');
            % purs_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, purs_spks, num2cell(Tpurs.pursuitOnset(Tpurs.angle==pur_best_dir)), 'uni', 1),'omitnan');
            % 
            % VMIs_best_eachEpo(c) = (stim_hz - sacc_hz) / (stim_hz + sacc_hz);
            % SPIs_best_eachEpo(c) = (sacc_hz - purs_hz) / (sacc_hz + purs_hz);
            % 
            % % each epoch separately
            % for d = 1:numel(dirs)
            %     mdir_spks = cellfun(@(q) q{c}, Tmdir.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tmdir.angle==dirs(d)), 'uni', 0);
            %     purs_spks = cellfun(@(q) q{c}, Tpurs.(sprintf('spiketimes_imec%d',Tclus.imec(c)))(Tpurs.angle==dirs(d)), 'uni', 0);
            % 
            %     stim_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=50 & (q-t(1))<150)/100)*1000, mdir_spks, Tmdir.TARG_ON(Tmdir.angle==dirs(d)), 'uni', 1),'omitnan');
            %     sacc_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, mdir_spks, num2cell(Tmdir.SACCADE(Tmdir.angle==dirs(d))), 'uni', 1),'omitnan');
            %     purs_hz = mean(cellfun(@(q,t) (sum((q-t(1))>=-50 & (q-t(1))<50)/100)*1000, purs_spks, num2cell(Tpurs.pursuitOnset(Tpurs.angle==dirs(d))), 'uni', 1),'omitnan');
            % 
            %     VMIs_eachDir{c,d} = (stim_hz - sacc_hz) / (stim_hz + sacc_hz);
            %     SPIs_eachDir{c,d} = (sacc_hz - purs_hz) / (sacc_hz + purs_hz);
            % end
        end
    
        Tclus.VMI_all = VMIs_all;
        Tclus.SPI_all = SPIs_all;
        Tclus.VMI_all_dp = VMIs_all_dp;
        Tclus.SPI_all_dp = SPIs_all_dp;
        Tclus.VMI_all_bc = VMIs_all_bc;
        Tclus.SPI_all_bc = SPIs_all_bc;
        % Tclus.VMI_best_highEpo = VMIs_best_highEpo;
        % Tclus.SPI_best_highEpo = SPIs_best_highEpo;
        % Tclus.VMI_best_eachEpo = VMIs_best_eachEpo;
        % Tclus.SPI_best_eachEpo = SPIs_best_eachEpo;
        % Tclus.VMI_eachDir = VMIs_eachDir;
        % Tclus.SPI_eachDir = SPIs_eachDir;

        % 
    
        if contains(string(Tclus.sessionName(1)),"scrappy")
            Tclus.monkey = repmat("scrappy",height(Tclus),1);
        else
            Tclus.monkey = repmat("yakko",height(Tclus),1);
        end
        Tclus.monkey = categorical(Tclus.monkey);

        CLUSTS_FILT = [CLUSTS_FILT; Tclus];
    end
end

% --- clusters --- %

% re-calculating SPI and VMI using only "good" trials 

%% FIGURE 5: SPI/VMI Distributions

f5a = figure;
f5a.Position = [100 100 1500 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

C = CLUSTS_FILT(~isnan(CLUSTS_FILT.VMI_all) & ~isnan(CLUSTS_FILT.SPI_all) & ~isinf(CLUSTS_FILT.VMI_all_dp) & ~isinf(CLUSTS_FILT.SPI_all_dp),:);
C = C(C.mdir1_Hz > FR_thresh & C.purs1_Hz > FR_thresh,:);
C = C((C.vis_sel_dir > seldir_range(1) & C.vis_sel_dir < seldir_range(2)) & (C.sac_sel_dir > seldir_range(1) & C.sac_sel_dir < seldir_range(2)) & (C.pur_sel_dir > seldir_range(1) & C.pur_sel_dir < seldir_range(2)),:);

ax1(1) = nexttile;
values = C.VMI_all_bc(C.monkey=='scrappy');
histStyle_KKN(values, 'BIN_WIDTH', 0.05, 'Y_LIMITS', [0 3600], 'X_LIMITS', [-1,1]);

ax1(2) = nexttile;
values = C.SPI_all_bc(C.monkey=='scrappy');
histStyle_KKN(values, 'BIN_WIDTH', 0.05, 'Y_LIMITS', [0 3600], 'X_LIMITS', [-1,1]);

ax2(1) = nexttile;
values = C.VMI_all_bc(C.monkey=='yakko');
histStyle_KKN(values, 'BIN_WIDTH', 0.05, 'Y_LIMITS', [0 1100], 'X_LIMITS', [-1,1]);

ax2(2) = nexttile;
values = C.SPI_all_bc(C.monkey=='yakko');
histStyle_KKN(values, 'BIN_WIDTH', 0.05, 'Y_LIMITS', [0 1100], 'X_LIMITS', [-1,1]);

%linkaxes(ax1, 'xy');
%linkaxes(ax2, 'xy');

%savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/spi_vmi_dists.pdf')

%%

f5b = figure;
f5b.Position = [100 100 750 800];
tl = tiledlayout(2,1);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

% Bin edges and centers
bin_edges = -1:0.1:1;
bin_centers = bin_edges(1:end-1) + 0.05;

% First monkey: SCRAPPY
ax(1) = nexttile;
monkey_name = 'scrappy';
C1 = C(C.monkey == monkey_name, :);

[counts, ~, ~] = histcounts2(C1.VMI_all, C1.SPI_all, bin_edges, bin_edges);
imagesc(bin_centers, bin_centers, counts');
set(gca, 'YDir', 'normal');
axis square;
colormap(gca, gray);  % white = more clusters
colorbar;
xlabel('VMI'); ylabel('SPI');
[r, p] = corr(C1.VMI_all, C1.SPI_all, 'rows','complete','Type','Pearson');
title(sprintf('%s: r = %0.3f, p = %0.3f', monkey_name, r, p));
prettyFig;

% Second monkey: YAKKO
ax(2) = nexttile;
monkey_name = 'yakko';
C2 = C(C.monkey == monkey_name, :);

[counts, ~, ~] = histcounts2(C2.VMI_all, C2.SPI_all, bin_edges, bin_edges);
imagesc(bin_centers, bin_centers, counts');
set(gca, 'YDir', 'normal');
axis square;
colormap(gca, gray);  % white = more clusters
colorbar;
xlabel('VMI'); ylabel('SPI');
[r, p] = corr(C2.VMI_all, C2.SPI_all, 'rows','complete','Type','Pearson');
title(sprintf('%s: r = %0.3f, p = %0.3f', monkey_name, r, p));
prettyFig;

savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/VMI_SPI_count_heatmap.pdf')

%% FIGURE 5:
f5b = figure;
f5b.Position = [100 100 750 800];
tl = tiledlayout(2,1);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

% First monkey
ax(1) = nexttile;
vmi_vals = C.VMI(C.monkey=='scrappy');
spi_vals = C.SPI(C.monkey=='scrappy');

plot(vmi_vals, spi_vals, 'ko');
hold on;
pfit = polyfit(vmi_vals, spi_vals, 1);
xfit = linspace(min(vmi_vals), max(vmi_vals), 100);
yfit = polyval(pfit, xfit);
plot(xfit, yfit, 'k-', 'LineWidth', 1.5);
[r,p] = corr(vmi_vals, spi_vals, 'Type', 'Pearson');
title(sprintf('Scrappy: r = %0.3f, p = %0.3f', r, p));
axis square
prettyFig;

% Second monkey
ax(2) = nexttile;
vmi_vals = C.VMI(C.monkey=='yakko');
spi_vals = C.SPI(C.monkey=='yakko');

plot(vmi_vals, spi_vals, 'ko');
hold on;
pfit = polyfit(vmi_vals, spi_vals, 1);
xfit = linspace(min(vmi_vals), max(vmi_vals), 100);
yfit = polyval(pfit, xfit);
plot(xfit, yfit, 'k-', 'LineWidth', 1.5);
[r,p] = corr(vmi_vals, spi_vals, 'Type', 'Pearson');
title(sprintf('Yakko: r = %0.3f, p = %0.3f', r, p));
axis square
prettyFig;

savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/VMI_v_SPI_scatter.pdf')

%% FIGURE 5: VMI and SPI as function of depth?

C1 = C(C.monkey == 'scrappy' & C.imec==0, :);
sessions = unique(C1.sessionName);
C_probe_depths = {[13,6];[7.2,9.8];[4.8,6.2];[8,11.5];[4.1,5.5];[11.5,9];[6.9,8.7];[7.7,9];[14,7.7];[10.7,8.4];[7.1,8.2];[8,8];[7.9,11.3];[8,11.3]};

VMI = []; SPI = []; depths = []; vis_sel = [];
for s = 1:numel(sessions)
    if ismember(s,[7]) %{9,10,11,1}, {1}, {2,12}, {3,7,8,13,14}
        VMI = [VMI; C1.VMI_all(C1.sessionName == sessions(s))];
        SPI = [SPI; C1.SPI_all(C1.sessionName == sessions(s))];
        depths = [depths; (C1.y_pos(C1.sessionName == sessions(s)) ./ 1000)]; % - C_probe_depths{s}(1)];
        vis_sel = [vis_sel; cos(deg2rad(C1.pur_pref_dir(C1.sessionName == sessions(s))))];
    end
end

% Set number of bins
n_bins = 10;

y_max = 6.2;
y_min = 0;

% Quantile-based bin edges
bin_edges = quantile(depths, linspace(0, 1, n_bins + 1));
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

mean_vmi_per_bin = NaN(size(bin_centers));
mean_spi_per_bin = NaN(size(bin_centers));
mean_vsel_per_bin = NaN(size(bin_centers));
cluster_counts = zeros(size(bin_centers));

for b = 1:length(bin_centers)
    in_bin = depths >= bin_edges(b) & depths < bin_edges(b + 1);
    if b == length(bin_centers)
        in_bin = depths >= bin_edges(b) & depths <= bin_edges(b + 1);
    end
    mean_vmi_per_bin(b) = mean(VMI(in_bin), 'omitnan');
    mean_spi_per_bin(b) = mean(SPI(in_bin), 'omitnan');
    mean_vsel_per_bin(b) = mean(vis_sel(in_bin), 'omitnan');
    cluster_counts(b) = sum(in_bin);
end

% --- Add empty bins at top and bottom ---
% Prepend and append to bin_edges
bin_edges = [y_min, bin_edges, y_max];
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% Pad metric vectors with NaN (empty bins)
mean_vmi_per_bin = [NaN, mean_vmi_per_bin, NaN];
mean_spi_per_bin = [NaN, mean_spi_per_bin, NaN];
mean_vsel_per_bin = [NaN, mean_vsel_per_bin, NaN];
cluster_counts    = [0, cluster_counts, 0];

% --- Create matrices for pcolor ---
x = [0 1];  % one column
vmi_matrix = repmat([mean_vmi_per_bin NaN]', 1, 2);
spi_matrix = repmat([mean_spi_per_bin NaN]', 1, 2);
vsel_matrix = repmat([mean_vsel_per_bin NaN]', 1, 2);

% --- Color maps ---
vmi_cmap = customColormap([133, 255, 19], [19, 133, 255], 256);
spi_cmap = customColormap([115, 30, 245], [245, 50, 30], 256);

% --- Plot VMI ---
f5c = figure;
f5c.Position = [100 100 600 1000];
subplot(1, 2, 1)
pcolor(x, bin_edges, vmi_matrix);
shading flat;
colormap(gca, vmi_cmap);
caxis([-0.5,0.5])
colorbar;
xticks([]);
yticks(bin_edges);
yticklabels(string(round(bin_edges, 1)));
ylim([y_min, y_max]);
prettyFig;

% --- Add cluster count labels (skip padding bins) ---
for b = 2:(n_bins+1)
    text(0.5, bin_centers(b), num2str(cluster_counts(b)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 10, 'Color', 'k');
end

% --- Plot SPI ---
subplot(1, 2, 2)
pcolor(x, bin_edges, spi_matrix);
shading flat;
colormap(gca, spi_cmap);
caxis([-0.5,0.5])
colorbar;
xticks([]);
yticks(bin_edges);
yticklabels(string(round(bin_edges, 1)));
ylim([y_min, y_max]);
prettyFig;

% --- Plot SPI ---
% subplot(1, 3, 3)
% pcolor(x, bin_edges, vsel_matrix);
% shading flat;
% colormap(gca, vmi_cmap);
% caxis([-1,1])
% colorbar;
% xticks([]);
% yticks(bin_edges);
% yticklabels(string(round(bin_edges, 1)));
% ylim([y_min, y_max]);
% prettyFig;

savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/VMI_SPI_depth_heatmap_142.pdf')



%% NEW THRESHOLDS FOR RSC
msOff_thresh  =  -50; % ms
sacLat_range =  [100 300];  % ms

sessions = fieldnames(PURS);

TBL = [];
for f = 1:numel(sessions)
    % --- mdir --- %
    Tmdir = MDIR.(sessions{f});
    Tmdir = Tmdir(Tmdir.result=='CORRECT',:);

    Tmdir.saccadeLatency = Tmdir.SACCADE-cellfun(@(q) q(1), Tmdir.FIX_OFF);
    Tmdir = Tmdir(Tmdir.saccadeLatency >= sacLat_range(1) & Tmdir.saccadeLatency < sacLat_range(2),:);
    
    [th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tmdir.eyePos, 'uni', 0);
    Tmdir = Tmdir(cellfun(@(r,t,s,a) (max(r(s-10:s+75)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s-10:s+75))) - a)) < 45), rh, th, num2cell(Tmdir.SACCADE), num2cell(Tmdir.angle), 'uni', 1),:);

    % --- purs --- %
    Tpurs = PURS.(sessions{f});
    Tpurs = Tpurs(Tpurs.result=='CORRECT',:);
    if contains(string(Tpurs.sessionName(1)), "Ya")
        Tpurs = Tpurs(isnan(Tpurs.csTimes(:,1) - Tpurs.pursuitOnset) | (Tpurs.csTimes(:,1) - Tpurs.pursuitOnset) < -50 | (Tpurs.csTimes(:,1) - Tpurs.pursuitOnset) > 50,:);
    else
        Tpurs = Tpurs(Tpurs.pursType=='pure',:);
        Tpurs = Tpurs(isnan(Tpurs.msOffset) | Tpurs.msOffset<msOff_thresh,:);

        [th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tpurs.eyePos, 'uni', 0);
        Tpurs = Tpurs(cellfun(@(r,t,s,a) (max(r(s+200:s+800)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s+200:s+800))) - a)) < 160), rh, th, num2cell(Tpurs.pursuitOnset), num2cell(Tpurs.angle), 'uni', 1),:);
    end
   
    new_row = table(Tmdir.sessionName(1), {Tmdir}, {Tpurs}, 'VariableNames', {'session','mdir','purs'});
    TBL = [TBL; new_row];
end

%% compiling all sessions together for rsc 
% TBL  =  mdir and purs for all sessions
% CLU  =  clusters across all sessions

FR_thresh     =  1;    % Hz
seldir_range =  [0.01 0.99]; 

CLU = CLUSTS_FILT(~isnan(CLUSTS_FILT.VMI_all) & ~isnan(CLUSTS_FILT.SPI_all) & ~isinf(CLUSTS_FILT.VMI_all_dp) & ~isinf(CLUSTS_FILT.SPI_all_dp),:);
CLU = CLU(CLU.mdir1_Hz > FR_thresh & CLU.purs1_Hz > FR_thresh,:);
CLU = CLU((CLU.vis_sel_dir > seldir_range(1) & CLU.vis_sel_dir < seldir_range(2)) & (CLU.sac_sel_dir > seldir_range(1) & CLU.sac_sel_dir < seldir_range(2)) & (CLU.pur_sel_dir > seldir_range(1) & CLU.pur_sel_dir < seldir_range(2)),:);
CLU = CLU(CLU.first_spike_sec<120,:);

sessions = unique(CLU.sessionName);
sessions = sessions([1,2,4,5,7:end]);

SPK_MATS = table();
for s = 1:length(sessions)
    fprintf('processing %s\n',sessions(s))
    this_sess = sessions(s);

    Tmdir = TBL.mdir{TBL.session == this_sess};
    Tpurs = TBL.purs{TBL.session == this_sess};

    imecs = unique(CLU.imec(CLU.sessionName == this_sess));
    for i = 1:numel(imecs)
        this_imec = imecs(i);
        Tclus = CLU(CLU.sessionName == this_sess & CLU.imec == this_imec,:);

        % mdir
        mdir_spks = cellfun(@(q,c) q(Tclus.cluster_id + 1), Tmdir.(sprintf('spiketimes_imec%d',this_imec)), 'uni', 0);
        sac_frs = cellfun(@(r,t) cellfun(@(q) (sum((q-t(1))>=-150 & (q-t(1))<=50) / 200) * 1000, r, 'uni', 1), mdir_spks, num2cell(Tmdir.SACCADE), 'uni', 0);
        sac_frs = vertcat(sac_frs{:}); % trials x neurons

        good_trls = mean(sac_frs,2) > (mean(mean(sac_frs,2)) - std(mean(sac_frs,2))*3) & mean(sac_frs,2) < (mean(mean(sac_frs,2)) + std(mean(sac_frs,2))*3);
        sac_frs = sac_frs(good_trls,:);
        Tmdir = Tmdir(good_trls,:);

        [~,~, mdir_conds] = unique([Tmdir.angle Tmdir.distance], 'rows');
        conds = unique(mdir_conds);

        z_sac_frs = NaN(size(sac_frs));  % Initialize the output
        for c = 1:numel(conds)
            cond_idx = mdir_conds == conds(c);       % Trials for this condition
            data_cond = sac_frs(cond_idx, :);        % Subset: [nTrials_in_cond x 238]
            
            % Z-score across trials (for each neuron = column)
            mu = mean(data_cond, 1, 'omitnan');      % Mean across trials
            sigma = std(data_cond, 0, 1, 'omitnan'); % Std across trials
            z_data_cond = (data_cond - mu) ./ sigma;
            
            % Assign back to the output
            z_sac_frs(cond_idx, :) = z_data_cond;
        end

        purs_spks = cellfun(@(q,c) q(Tclus.cluster_id + 1), Tpurs.(sprintf('spiketimes_imec%d',this_imec)), 'uni', 0);
        pur_frs = cellfun(@(r,t) cellfun(@(q) (sum((q-t(1))>=-50 & (q-t(1))<=150) / 200) * 1000, r, 'uni', 1), purs_spks, num2cell(Tpurs.pursuitOnset), 'uni', 0);
        pur_frs = vertcat(pur_frs{:}); % trials x neurons

        good_trls = mean(pur_frs,2) > (mean(mean(pur_frs,2)) - std(mean(pur_frs,2))*3) & mean(pur_frs,2) < (mean(mean(pur_frs,2)) + std(mean(pur_frs,2))*3);
        pur_frs = pur_frs(good_trls,:);
        Tpurs = Tpurs(good_trls,:);

        [~,~, purs_conds] = unique([Tpurs.angle Tpurs.pursuitSpeed], 'rows');
        conds = unique(purs_conds);

        z_pur_frs = NaN(size(pur_frs));  % Initialize the output
        for c = 1:numel(conds)
            cond_idx = purs_conds == conds(c);       % Trials for this condition
            data_cond = pur_frs(cond_idx, :);        % Subset: [nTrials_in_cond x 238]
            
            % Z-score across trials (for each neuron = column)
            mu = mean(data_cond, 1, 'omitnan');      % Mean across trials
            sigma = std(data_cond, 0, 1, 'omitnan'); % Std across trials
            z_data_cond = (data_cond - mu) ./ sigma;
            
            % Assign back to the output
            z_pur_frs(cond_idx, :) = z_data_cond;
        end

        if contains(string(Tclus.sessionName(1)),"scrappy")
            monkey = categorical("scrappy");
        else
            monkey = categorical("yakko");
        end

        good_clusts = sum(isnan(z_sac_frs),1)==0 & sum(isnan(z_pur_frs),1)==0;

        new_row = table(monkey, this_sess, this_imec, {Tclus(good_clusts',:)}, {z_sac_frs(:,good_clusts)}, {z_pur_frs(:,good_clusts)}, 'VariableNames', {'monkey','session','imec','clusts','mdir','purs'});
        SPK_MATS = [SPK_MATS; new_row];
    end
end


%% calculating rsc values for all pairs

RSCs = cell(sum(cellfun(@(q) height(q) * (height(q) - 1) / 2, SPK_MATS.clusts, 'uni', 1))*2, 15);
counter = 1;
for s = 1:height(SPK_MATS)
    fprintf('Processing session: %s\n', SPK_MATS.session(s));

    C = SPK_MATS.clusts{s};
    Tmdir = SPK_MATS.mdir{s};
    Tpurs = SPK_MATS.purs{s};

    coords = [C.x_pos C.y_pos];

    for i = 1:height(C)%-1
        for j = 1:height(C)%i+1:height(C)
            n_dist = sqrt((coords(i,1) - coords(j,1))^2 + (coords(i,2) - coords(j,2))^2);   
            if n_dist >= 50
                % mdir 
                sac_dir_diff = min(abs(C.sac_pref_dir(i) - C.sac_pref_dir(j)), 360 - abs(C.sac_pref_dir(i) - C.sac_pref_dir(j)));
                [rho_sac,pval_sac] = corr(Tmdir(:,i), Tmdir(:,j));

                % purs
                pur_dir_diff = min(abs(C.pur_pref_dir(i) - C.pur_pref_dir(j)), 360 - abs(C.pur_pref_dir(i) - C.pur_pref_dir(j)));
                [rho_pur,pval_pur] = corr(Tpurs(:,i), Tpurs(:,j));

                RSCs(counter, :) = {SPK_MATS.monkey(s), sprintf('%s_%d',string(SPK_MATS.session(s)),SPK_MATS.imec(s)), C.cluster_id(i), ...
                                    C.cluster_id(j), C.VMI_all(i), C.SPI_all(i), C.VMI_all(j), C.SPI_all(j), ...
                                    n_dist, sac_dir_diff, pur_dir_diff, rho_sac, pval_sac, rho_pur, pval_pur};

                counter = counter + 1;
            end    
        end
    end 
end

RSCs(any(cellfun(@isempty, RSCs), 2), :) = [];
Rtbl = cell2table(RSCs, 'VariableNames', {'monkey', 'session', 'c1_id', 'c2_id', 'c1_vmi', 'c1_spi', 'c2_vmi', 'c2_spi', 'dist', 'sac_dir_diff', 'pur_dir_diff', 'rho_sac', 'pval_sac', 'rho_pur', 'pval_pur'});
Rtbl.session = categorical(Rtbl.session);

%% FIGURE 6: Histograms of rsc values 

f6a = figure;
f6a.Position = [100 100 1800 500];
tl = tiledlayout(1,4);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';


ax1(1) = nexttile;
values = rtoZ(Rtbl.rho_sac(Rtbl.monkey=='scrappy'));
histStyle_KKN(values, 'BIN_WIDTH', 0.1); %, 'Y_LIMITS', [0 3600], 'X_LIMITS', [-1,1]);

ax2(1) = nexttile;
values = rtoZ(Rtbl.rho_sac(Rtbl.monkey=='yakko'));
histStyle_KKN(values, 'BIN_WIDTH', 0.1); %, 'Y_LIMITS', [0 1100], 'X_LIMITS', [-1,1]);

ax1(2) = nexttile;
values = rtoZ(Rtbl.rho_pur(Rtbl.monkey=='scrappy'));
histStyle_KKN(values, 'BIN_WIDTH', 0.1); %, 'Y_LIMITS', [0 3600], 'X_LIMITS', [-1,1]);

ax2(2) = nexttile;
values = rtoZ(Rtbl.rho_pur(Rtbl.monkey=='yakko'));
histStyle_KKN(values, 'BIN_WIDTH', 0.1); %, 'Y_LIMITS', [0 1100], 'X_LIMITS', [-1,1]);


%% FIGURE 6: rsc as fxn of distance between contacts 

N_BINS = 7;

% even width 
bin_edges = linspace(50,4500,N_BINS);
bin_edges = bin_edges./1000;

bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

f6a = figure;
f6a.Position = [100 100 1800 500];
tl = tiledlayout(1,4);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

% scrappy
sessions = unique(Rtbl.session(Rtbl.monkey == 'scrappy'));
sessions = sessions([6,10,14,15,17,21,23,24,25,26]);

ax1(1) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'scrappy' & Rtbl.session == sessions(s), :);
    Rmonk.dist = Rmonk.dist ./ 1000;  % Convert to mm

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.dist >= bin_edges(b) & Rmonk.dist < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_sac(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 4.5])
    ylim([-0.018, 0.062])
    hold on;
end
prettyFig;
axis square

% Combine data across all selected sessions
Rscr_sac = Rtbl(Rtbl.monkey == 'scrappy' & ismember(Rtbl.session, sessions), :);
Rscr_sac.dist = Rscr_sac.dist ./ 1000;
valid_rows = ~isnan(Rscr_sac.rho_sac);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_sac.dist(valid_rows), rtoZ(Rscr_sac.rho_sac(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(1));

ax1(3) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'scrappy' & Rtbl.session == sessions(s), :);
    Rmonk.dist = Rmonk.dist ./ 1000;  % Convert to mm

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.dist >= bin_edges(b) & Rmonk.dist < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_pur(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 4.5])
    ylim([-0.018, 0.062])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_pur = Rtbl(Rtbl.monkey == 'scrappy' & ismember(Rtbl.session, sessions), :);
Rscr_pur.dist = Rscr_pur.dist ./ 1000;
valid_rows = ~isnan(Rscr_pur.rho_pur);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_pur.dist(valid_rows), rtoZ(Rscr_pur.rho_pur(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(3));

sessions = unique(Rtbl.session(Rtbl.monkey == 'yakko'));

ax1(2) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'yakko' & Rtbl.session == sessions(s), :);
    Rmonk.dist = Rmonk.dist ./ 1000;  % Convert to mm

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.dist >= bin_edges(b) & Rmonk.dist < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_sac(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 4.5])
    ylim([-0.018, 0.062])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_sac = Rtbl(Rtbl.monkey == 'yakko', :);
Rscr_sac.dist = Rscr_sac.dist ./ 1000;
valid_rows = ~isnan(Rscr_sac.rho_sac);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_sac.dist(valid_rows), rtoZ(Rscr_sac.rho_sac(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(2));

ax1(4) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'yakko' & Rtbl.session == sessions(s), :);
    Rmonk.dist = Rmonk.dist ./ 1000;  % Convert to mm

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.dist >= bin_edges(b) & Rmonk.dist < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_pur(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 4.5])
    ylim([-0.018, 0.062])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_pur = Rtbl(Rtbl.monkey == 'yakko', :);
Rscr_pur.dist = Rscr_pur.dist ./ 1000;
valid_rows = ~isnan(Rscr_pur.rho_pur);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_pur.dist(valid_rows), rtoZ(Rscr_pur.rho_pur(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(4));

savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/rsc_dist.pdf')


%% FIGURE 6: rsc as fxn of diff in pref direction

N_BINS = 7;

% even width 
bin_edges = linspace(0,180,N_BINS);

bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

f6a = figure;
f6a.Position = [100 100 1800 500];
tl = tiledlayout(1,4);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

% scrappy
sessions = unique(Rtbl.session(Rtbl.monkey == 'scrappy'));
sessions = sessions([6,10,14,15,17,21,23,24,25,26]);

ax1(1) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'scrappy' & Rtbl.session == sessions(s), :);

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.sac_dir_diff >= bin_edges(b) & Rmonk.sac_dir_diff < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_sac(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 190])
    ylim([-0.0075, 0.04])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_sac = Rtbl(Rtbl.monkey == 'scrappy' & ismember(Rtbl.session, sessions), :);
valid_rows = ~isnan(Rscr_sac.rho_sac);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_sac.sac_dir_diff(valid_rows), rtoZ(Rscr_sac.rho_sac(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(1));

ax1(3) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'scrappy' & Rtbl.session == sessions(s), :);

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.pur_dir_diff >= bin_edges(b) & Rmonk.pur_dir_diff < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_pur(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 190])
    ylim([-0.0075, 0.04])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_pur = Rtbl(Rtbl.monkey == 'scrappy' & ismember(Rtbl.session, sessions), :);
valid_rows = ~isnan(Rscr_pur.rho_pur);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_pur.pur_dir_diff(valid_rows), rtoZ(Rscr_pur.rho_pur(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(3));

sessions = unique(Rtbl.session(Rtbl.monkey == 'yakko'));

ax1(2) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'yakko' & Rtbl.session == sessions(s), :);

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.sac_dir_diff >= bin_edges(b) & Rmonk.sac_dir_diff < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_sac(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 190])
    ylim([-0.0075, 0.04])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_sac = Rtbl(Rtbl.monkey == 'yakko', :);
valid_rows = ~isnan(Rscr_sac.rho_sac);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_sac.sac_dir_diff(valid_rows), rtoZ(Rscr_sac.rho_sac(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(2));

ax1(4) = nexttile;
for s = 1:length(sessions)
    Rmonk = Rtbl(Rtbl.monkey == 'yakko' & Rtbl.session == sessions(s), :);

    mean_rho = NaN(1, N_BINS - 1);
    sem_rho = NaN(1, N_BINS - 1);

    for b = 1:length(bin_centers)
        in_bin = Rmonk.pur_dir_diff >= bin_edges(b) & Rmonk.pur_dir_diff < bin_edges(b + 1);
        n_pairs = sum(in_bin);
        if n_pairs < 100
            continue;  % Leave as NaN to skip but still connect line
        end

        rho_vals = rtoZ(Rmonk.rho_pur(in_bin));
        mean_rho(b) = mean(rho_vals, 'omitnan');
        sem_rho(b) = std(rho_vals, 'omitnan') / sqrt(n_pairs);
    end

    % Plot with NaNs still in place – MATLAB will connect valid points
    errorbar(bin_centers, mean_rho, sem_rho, 'ko-', 'LineWidth', 1.5);
    hold on;
    yline(0,'k--')
    xlim([0 190])
    ylim([-0.0075, 0.04])
    hold on;
end
axis square
prettyFig;

% Combine data across all selected sessions
Rscr_pur = Rtbl(Rtbl.monkey == 'yakko', :);
valid_rows = ~isnan(Rscr_pur.rho_pur);
n_pairs = sum(valid_rows);
[rval, pval] = corr(Rscr_pur.pur_dir_diff(valid_rows), rtoZ(Rscr_pur.rho_pur(valid_rows)), 'type', 'Pearson');

% Add text box
text(0.2, 0.05, sprintf('n = %d\nr = %.3f, p = %.3g', n_pairs, rval, pval), ...
    'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Parent', ax1(4));

savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/rsc_prefdirdiff.pdf')

%% 

f6a = figure;
f6a.Position = [100 100 1800 500];
tl = tiledlayout(1,4);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

edges = [-1:0.25:1];
bin_centers = edges(1:end-1) + diff(edges)/2;

%Rtbl2 = Rtbl;
%Rtbl2.c1_id = Rtbl.c2_id;
%Rtbl2.c2_id = Rtbl.c1_id;
%Rtbl3 = [Rtbl; Rtbl2];
Rtbl3 = Rtbl;

% Preallocate matrix for binned rho_sac values
heatmap_vals = NaN(length(edges)-1, length(edges)-1);

% scrappy
sessions = unique(Rtbl3.session(Rtbl3.monkey == 'scrappy'));
sessions = sessions([6,10,14,15,17,21,23,24,25,26]);

Rmonk = Rtbl3(Rtbl3.monkey == 'scrappy', :);

ax1(1) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_spi >= edges(i) & Rmonk.c1_spi < edges(i+1);
        in_y_bin = Rmonk.c2_spi >= edges(j) & Rmonk.c2_spi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_sac(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end



% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

ax1(3) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_spi >= edges(i) & Rmonk.c1_spi < edges(i+1);
        in_y_bin = Rmonk.c2_spi >= edges(j) & Rmonk.c2_spi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_pur(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

% scrappy
Rmonk = Rtbl3(Rtbl3.monkey == 'yakko', :);

ax1(2) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_spi >= edges(i) & Rmonk.c1_spi < edges(i+1);
        in_y_bin = Rmonk.c2_spi >= edges(j) & Rmonk.c2_spi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_sac(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

ax1(4) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_spi >= edges(i) & Rmonk.c1_spi < edges(i+1);
        in_y_bin = Rmonk.c2_spi >= edges(j) & Rmonk.c2_spi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_pur(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

%savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/rsc_spi_heatmap.pdf')

%%

f6a = figure;
f6a.Position = [100 100 1800 500];
tl = tiledlayout(1,4);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

edges = [-1:0.25:1];
bin_centers = edges(1:end-1) + diff(edges)/2;

% Preallocate matrix for binned rho_sac values
heatmap_vals = NaN(length(edges)-1, length(edges)-1);

% scrappy
sessions = unique(Rtbl.session(Rtbl.monkey == 'scrappy'));
%sessions = sessions([6,10,14,15,17,21,23,24,25,26]);

sessions = sessions([1,3,8,11,12,13,14,15,17,18,19,20,21]); %1,3,8,11,12,13, 14*,15,17,18,19*,20,21,

Rmonk = Rtbl3(Rtbl3.monkey == 'scrappy' & ismember(Rtbl3.session,sessions), :);

ax1(1) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_vmi >= edges(i) & Rmonk.c1_vmi < edges(i+1);
        in_y_bin = Rmonk.c2_vmi >= edges(j) & Rmonk.c2_vmi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_sac(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

ax1(3) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_vmi >= edges(i) & Rmonk.c1_vmi < edges(i+1);
        in_y_bin = Rmonk.c2_vmi >= edges(j) & Rmonk.c2_vmi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_pur(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

% scrappy
Rmonk = Rtbl3(Rtbl3.monkey == 'yakko', :);

ax1(2) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_vmi >= edges(i) & Rmonk.c1_vmi < edges(i+1);
        in_y_bin = Rmonk.c2_vmi >= edges(j) & Rmonk.c2_vmi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_sac(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

ax1(4) = nexttile;
% Loop over bins
for i = 1:length(edges)-1  % x-axis: c1_spi
    for j = 1:length(edges)-1  % y-axis: c2_spi
        in_x_bin = Rmonk.c1_vmi >= edges(i) & Rmonk.c1_vmi < edges(i+1);
        in_y_bin = Rmonk.c2_vmi >= edges(j) & Rmonk.c2_vmi < edges(j+1);
        in_bin = in_x_bin & in_y_bin;
        
        if any(in_bin)
            heatmap_vals(j, i) = mean(Rmonk.rho_pur(in_bin), 'omitnan');  % note: row = y, col = x
        end
    end
end

% Plot with pcolor
pcolor(edges, edges, padarray(heatmap_vals, [1 1], NaN, 'post'));  % pad for pcolor alignment
shading flat;
colormap(parula);  % choose color map you like
caxis([-0.06 0.06]);
colorbar;
axis square;
prettyFig;

%savebigPDF(1, '/Users/kendranoneman/Posters/GRC-2025/rsc_vmi_heatmap.pdf')

