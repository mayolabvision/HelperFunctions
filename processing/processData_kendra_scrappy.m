clear
clc

addpath(genpath('/Users/kendranoneman/Packages')) % add nevUtils and HelperFunctions to path
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paths 
NET_PATH   =  '/Users/kendranoneman/Packages/nasnet/networks';
DATA_PATH  =  '/Users/kendranoneman/OneDrive/DATA';
CSV_PATH   =  '/Users/kendranoneman/OneDrive/DATA/RECORDING_INFO.csv';

% Session details
EXPERIMENTER  =  'kendra';
MONKEY        =  'scrappy';
SESSION       =  '0097a';

% Spike thresholding
GAMMA  =  0.2; 

% RF Mapping
INTERP_RF = false; 
FIRST_BIN = 0; BIN_WIDTH = 50; BIN_STEP = 10; NBINS = 24;

% Smooth Pursuit
PURS_PREINT  =  25;
PURS_POSTINT =  210;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load in recording details and make output directory
[this_sess,filename]  =  read_recordingNotes(CSV_PATH,EXPERIMENTER,MONKEY,SESSION);

taskTypes = {'rfmp','purs','mdir','fstm'};
numTasks = [this_sess.rfmp_num this_sess.purs_num this_sess.mdir_num this_sess.fstm_num];
tasks = [];
for i=1:length(taskTypes)
    tasks = [tasks, arrayfun(@(x) sprintf('%s%d', taskTypes{i}, x), 1:numTasks(i), 'UniformOutput', false)];
end

output_path = fullfile(DATA_PATH, 'processed', filename);
if ~exist(output_path, 'dir'), mkdir(output_path); end
if ~exist(fullfile(output_path, 'figs'), 'dir'), mkdir(fullfile(output_path, 'figs')); end

%% Extracting raw data from nev/out datafiles 
if ~exist(fullfile(output_path, 'dataTable666.mat'), 'file')
    % Name/organize channels based on recording details
    [mappings,probe_specs] = map_channelsNumbersToNames(this_sess.mapFile_name,this_sess.probeID{1},'probeDepths_mm',this_sess.recordDepth_mm{1});

    tic

    % Make structure to hold all data 
    S1 = struct();
    S1.recording_info = table2struct(this_sess);
    S1.channels = mappings;
    
    for i = 1:length(taskTypes)
        these_tasks = arrayfun(@(x) sprintf('%s%d', taskTypes{i}, x), 1:numTasks(i), 'UniformOutput', false);
        for f = 1:length(these_tasks)
            this_task = these_tasks{f};

            % Step 1. Pull out data, spike sorting
            fprintf('\n---- generating nev_out for %s ----\n', this_task);
        
            nevname = sprintf('%s/raw/%s_%s.nev', DATA_PATH, filename, this_task);
            [nev, out] = extract_nevout(nevname, 'SPIKE_SORT', true, 'netFolder', NET_PATH);
        
            % Step 2. Extract raw data and waveforms
            % if ~exist(fullfile(output_path, sprintf('%s_rawData.mat',this_task)), 'file')
            %     fprintf('\n---- extracting raw/waveforms for %s ----\n', this_task);
            %     raw = extract_rawData(nev,out);
            %     raw = cellfun(@(q) single(q), raw, 'uni', 0); 
            %     save(fullfile(output_path, sprintf('%s_rawData.mat',this_task)), 'raw', '-v7.3');
            %     clear raw;
            % end

            % Step 3. Generate and save data table to structure
            fprintf('\n---- generating table for %s ----\n', this_task);
            tbl = format_dataTable(nev, out, this_task);

            if ismember('IGNORED', tbl.Properties.VariableNames)
                tbl = movevars(tbl,{'IGNORED'},'After','CORRECT');
            end

            % Convert structures to a cell array of string representations
            all_params = {tbl.params.block}.';
            structStrings = cellfun(@(x) jsonencode(x), all_params, 'UniformOutput', false);
            [~, uniqueIdx] = unique(structStrings, 'stable');
            unique_structs = all_params(uniqueIdx);
            merged_struct = struct();
            fieldNames = fieldnames(unique_structs{1});
            for ii = 1:numel(fieldNames)
                field = fieldNames{ii};
                merged_struct.(field) = cellfun(@(s) s.(field), unique_structs, 'UniformOutput', false);
                if all(cellfun(@isnumeric, merged_struct.(field)))
                    merged_struct.(field) = cell2mat(merged_struct.(field));
                end
            end
            tbl.params = [];

            S1.(this_task).hdr = out.hdr;
            S1.(this_task).params = merged_struct;
            S1.(this_task).data = tbl;
            clear tbl;
        end
    end

    S = unify_taskTables(S1,taskTypes);

    % Save the structure S to the specified file
    save(fullfile(output_path, 'dataTable.mat'), 'S');
    
    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

else
    fprintf('\n---- loading S for %s ----\n', filename);
    load(fullfile(output_path, 'dataTable.mat'))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PARAMETERS

num_channels  =  height(S.channels);
good_chans    =  S.channels.ripChan_num(S.channels.depth_order>0);
num_rfmp = this_sess.rfmp_num; num_purs = this_sess.purs_num; num_mdir = this_sess.mdir_num; num_fstm = this_sess.fstm_num;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TEST

T = S.mdir1.data;

f = figure;
plot(T.eyePos(:,1),'k-')
hold on;
xline(T.FIX_ON)


%% ------------------------------------------- 1. Overall Performance -------------------------------------------
fig_path = fullfile(output_path, 'figs', 'session_trialOutcomes.png');

main_taskTypes = taskTypes(cellfun(@(q) ~contains(q,'fstm'), taskTypes, 'uni', 1));
main_tasks = tasks(cellfun(@(q) ~contains(q,'fstm'), tasks, 'uni', 1));
results = cell(length(main_tasks),1);
for task=1:length(main_tasks)
    results{task} = S.(main_tasks{task}).data.result;
end

f1a = figure;
f1a.Position = [100 100 1800 400];
tl = tiledlayout(1,length(main_tasks));
tl.TileSpacing = 'tight';
tl.Padding = 'compact';

set(gcf,'color','w')
for ii = 1:length(main_tasks)
    nexttile
    pieChart_trialOutcomes(results{ii}, main_tasks{ii}, vertcat(results{:}));
end
title(tl,sprintf('%s',filename),'fontsize',18,'interpreter','none')
print(f1a, fig_path, '-dpng', '-r300');

%% ------------------------------------------- 2. RF Mapping -------------------------------------------
fig_path = fullfile(output_path, 'figs', 'rfmp');
if ~exist(fullfile(output_path, 'figs', 'rfmp'), 'dir'), mkdir(fullfile(output_path, 'figs', 'rfmp')); end
rfmp_tasks = tasks(cellfun(@(q) contains(q,'rfmp'), tasks, 'uni', 1));


all_FRs = cell(num_channels,num_rfmp);
for batch = 1:num_rfmp
    T = S.(rfmp_tasks{batch}).data;

    all_FRs(:,batch) = format_tableToRFMap(T,FIRST_BIN,BIN_WIDTH,BIN_STEP,NBINS,GAMMA);
end

% RF Map plot for each unit
for batch = 1:num_rfmp
    for unit = 1:length(good_chans)
        f2a = figure('Visible','off');
        f2a.Position = [100 100 1800 900];

        chan_name = S.channels.mapped_name{good_chans(unit)};
        tl = heatMap_rfOverTime(all_FRs{good_chans(unit),batch},'BIN_EDGES',bin_edges, 'INTERP', INTERP_RF,...
                               'X_VALS',pix2deg(xvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)), ...
                               'Y_VALS',pix2deg(yvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)));
        
        title(tl,sprintf('%s_%s',filename,rfmp_tasks{batch}),'fontsize',20,'interpreter','none')
        subtitle(tl,sprintf('%s (ripChan = %d)',chan_name, S.channels.ripChan_num(good_chans(unit))),'fontsize',16,'interpreter','none')
    
        print(f2a, fullfile(fig_path, sprintf('%s-%s.png', rfmp_tasks{batch}, chan_name)), '-dpng', '-r200');
        fprintf(sprintf('\n----Unit %.2d complete----',unit))
    end
end
fprintf('\n----------------------\n')

% TODO: Compare RF maps, Quantify how different they are 

%% ------------------------------------------- 3. Smooth Pursuit -------------------------------------------
fig_path = fullfile(output_path, 'figs', 'purs');
if ~exist(fullfile(output_path, 'figs', 'purs'), 'dir'), mkdir(fullfile(output_path, 'figs', 'purs')); end

purs_tasks = tasks(cellfun(@(q) contains(q,'purs'), tasks, 'uni', 1));

T = [];
for batch = 1:num_purs
    T = [T; S.(purs_tasks{batch}).data];
end

% A. Checks for biases in percent correct across conditions
f3a = figure('Visible','off');
f3a.Position = [100 100 950 400];

tl = polarPlot_pursDistByCondition(T,'perCorrect');
title(tl,sprintf('%s_purs',filename),'fontsize',20,'interpreter','none')
subtitle(tl, {'';'% Correct'},'fontsize',14)

print(f3a, fullfile(fig_path, 'perCorrByCond.png'), '-dpng', '-r300');

% B. Percent pure pursuit trials across conditions
f3b = figure('Visible','off');
f3b.Position = [100 100 950 400];

tl = polarPlot_pursDistByCondition(T,'perPure');
title(tl,sprintf('%s_purs',filename),'fontsize',20,'interpreter','none')
subtitle(tl, {'';'% of Trials w/ Pure Pursuit Initiation'},'fontsize',14)

print(f3b, fullfile(fig_path, 'perPureByCond.png'), '-dpng', '-r300');

% C. Pursuit latencies across conditions
f3c = figure('Visible','off');
f3c.Position = [100 100 950 400];

tl = polarPlot_pursDistByCondition(T,'pursLatency');
title(tl,sprintf('%s_purs',filename),'fontsize',20,'interpreter','none')
subtitle(tl, {'';'Pursuit Latency (ms)'},'fontsize',14)

print(f3c, fullfile(fig_path, 'latencyByCond.png'), '-dpng', '-r300');


% Plot "pure pursuit" eye traces, split by speeds/jumps
rt = S.(purs_tasks{1}).params.reactionTime;

% Pure pursuit only
f3d = figure('Visible','off');
f3d.Position = [100 100 1500 900];

[tl,num_pure] = eyeTraces_pursSplitByConditions(T,rt,'PREINT',PURS_PREINT,'POSTINT',PURS_POSTINT,'PURE_ONLY',true);

title(tl,sprintf('%s_purs',filename),'fontsize',20,'interpreter','none')
subtitle(tl, sprintf('(# of pure pursuit trials = %d, step-ramp duration = %d ms)',num_pure,rt))

print(f3d, fullfile(fig_path, 'eyeTraces_pureTrials.png'), '-dpng', '-r300');

% All trials 
f3e = figure('Visible','off');
f3e.Position = [100 100 1500 900];

[tl,num_pure] = eyeTraces_pursSplitByConditions(T,rt,'PREINT',PURS_PREINT,'POSTINT',PURS_POSTINT,'PURE_ONLY',false);

title(tl,sprintf('%s_purs',filename),'fontsize',20,'interpreter','none')
subtitle(tl, sprintf('(# of pure pursuit trials = %d, step-ramp duration = %d ms)',num_pure,rt))

print(f3e, fullfile(fig_path, 'eyeTraces_allTrials.png'), '-dpng', '-r300');

%% Pursuit rasters

speed = 10; jump = -1; 
csFlag = 0;

unit = 16;

f3f = figure; %('Visible','off');
f3f.Position = [100 100 1800 900];

angles = sort(unique(T.angle))';
angle_order = [6,3,2,1,4,7,8,9];
y_lims = []; % Store y-axis limits
frs_perAng = cell(length(angles),1);
for ang = 1:length(angles)
    these_trls = T(T.angle==angles(ang) & T.pursuitSpeed==speed & T.jump==jump & T.csFlag==csFlag,:);

    sptimes = cellfun(@(q,z,r) q(z)-r, these_trls.spiketimes(:,unit), cellfun(@(q) q>GAMMA, these_trls.net_labels(:,unit), 'uni', 0), num2cell(these_trls.PURSUIT_TARG), 'uni', 0);

    subplot(3,3,angle_order(ang))

    line_color = [0,0,0]./255; sem_shade = [200,200,200]./255;
    raster_sdf(sptimes', [-25, 500], 10, 'line_color', line_color, 'sem_shade', sem_shade)
    ax = gca;
    y_lims = [y_lims; ax.YLim];
    frs_perAng{ang} = cellfun(@(q) (sum(q>=0 & q <500)*(1000/500)), sptimes, 'uni', 1);

    % Store y-axis limits
    

    blah = 1;
end

% Pad each cell with NaNs to match maxLength
maxLength = max(cellfun(@numel, frs_perAng));
frs_perAng = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_perAng, 'UniformOutput', false);

stimrate = vertcat(frs_perAng{:})';

% Generate randomized index of stimrate values, WITH REPLACEMENT
shuffles = 1000;
rhoPst = [];

for sh=1:shuffles
    randind=randi( (size(stimrate,1)*size(stimrate,2)), size(stimrate,1), size(stimrate,2) );
    permutedStimrate = stimrate(randind);
    rhoPst = [rhoPst; nanmean(permutedStimrate)];
end

sorted_rhoPst=sort(rhoPst);
rhoLst = sorted_rhoPst(shuffles*.05,:); % 95% lower confidence interval
rhoUst = sorted_rhoPst(shuffles-(shuffles*.05),:); % 95% upper confidence interval

% calculate tuning preferences
%theta = 0:360/length(a.CND):360; theta(end)=[];
theta = 0:45:315;
[visds, visdp] = tuningbias(theta,nanmean(stimrate));

subplot(3,3,5)
rho = nanmean(stimrate);
dst = sprintf('%0.2f',visds);
dpt = sprintf('%0.2f',visdp);
polarplot(deg2rad([theta 0]),[rho rho(1)],'ko-',...
    'markerfacecolor','k','linewidth',3)
hold on
polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'go--','LineWidth',2);
polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'go--','LineWidth',2);   

h1=polarplot(deg2rad(visdp),max(rho),'k^'); % Plot black triangle at best dir
txd1 = get(h1,'ThetaData');
tyd1 = get(h1,'RData');
set(h1,'ThetaData',txd1(1),'RData',tyd1(1));
set(h1,'markersize',10,'markerfacecolor','k'); hold off;
title(['VisDir: ',dpt,', Sel: ',dst]);
prettyFig


% Find the global y-axis limits
global_y_lim = [min(y_lims(:,1)), max(y_lims(:,2))];

% Apply the limits to all subplots
for ang = 1:length(angle_order)
    subplot(3,3,ang)
    if ang ~= 5
        ylim(global_y_lim);
    end
end

han=axes(f3f,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
xlabel(han,{'';'time aligned to target motion onset (ms)'},'fontsize',16);
title(han,{sprintf('%s_%s',filename,rfmp_tasks{batch}); sprintf('%s (ripChan = %d)',chan_name, S.channels.ripChan_num(good_chans(unit)))},'fontsize',20,'interpreter','none')
    
% print(f3f, fullfile(fig_path, sprintf('%s-%s.png', rfmp_tasks{batch}, chan_name)), '-dpng', '-r200');
fprintf(sprintf('\n----Unit %.2d complete----',unit))




