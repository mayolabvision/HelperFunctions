clear
clc

addpath(genpath('/Users/kendranoneman/Packages')) % add nevUtils and HelperFunctions to path
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions'))

%% SESSION PARAMETERS

netFolder = '/Users/kendranoneman/Packages/nasnet/networks';
dataFolder = '/Users/kendranoneman/OneDrive/DATA';
csvPath = '/Users/kendranoneman/OneDrive/DATA/RECORDING_INFO.csv';

EXPERIMENTER = 'kendra';
MONKEY = 'scrappy';
SESSION = '0077a';

% Load in details from notes
[this_sess,filename] = read_recordingNotes(csvPath,EXPERIMENTER,MONKEY,SESSION);

% Name/organize channels based on recording details
[mappings,probe_specs] = map_channelsNumbersToNames(this_sess.mapFile_name,this_sess.probeID{1},'probeDepths_mm',this_sess.recordDepth_mm{1});

%% Extracting/saving nev and out, using Smith lab nevutils

tic
taskTypes = {'rfmp','purs','mdir','fstm'};
numTasks = [this_sess.rfmp_num this_sess.purs_num this_sess.mdir_num this_sess.fstm_num];

% Make structure to hold all data 
S1 = struct();
S1.recording_info = table2struct(this_sess);
S1.channels = mappings;

for i = 1:length(taskTypes)
    if ~exist(fullfile(dataFolder, 'processed', sprintf('%s.mat', filename)), 'file')
        these_tasks = arrayfun(@(x) sprintf('%s%d', taskTypes{i}, x), 1:numTasks(i), 'UniformOutput', false);
        for f = 1:length(these_tasks)
            this_task = these_tasks{f};
            if ~exist(fullfile(dataFolder,'nev_out',sprintf('%s_%s.mat',filename,this_task)), 'file')
                fprintf('\n---- generating nev_out for %s ----\n', this_task);
            
                nevname = sprintf('%s/raw/%s_%s.nev', dataFolder, filename, this_task);
                [nev, out] = extract_nevout(nevname, 'SPIKE_SORT', true, 'netFolder', netFolder);
            
                save(fullfile(dataFolder,'nev_out',sprintf('%s_%s.mat',filename,this_task)),'nev','out')
            else
                fprintf('\n---- loading nev_out for %s ----\n', this_task);
                load(fullfile(dataFolder,'nev_out',sprintf('%s_%s.mat',filename,this_task)))
            end
        
            % Store the table data for this task
            fprintf('\n---- generating table for %s ----\n', this_task);
            [dat,tbl] = format_dataTable(nev, out, this_task);
            save(fullfile(dataFolder, 'dats', sprintf('%s_%s.mat', filename, this_task)), 'dat');
    
            S1.(this_task) = tbl;
        end
    end
end

S = unify_taskTables(S1,taskTypes);

% Save the structure S to the specified file
filePath = fullfile(dataFolder, 'processed', sprintf('%s.mat', filename));
save(filePath, 'S');

tc = toc;
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf(sprintf('Elapsed time is %2.2f minutes',tc/60))
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');


%% Fixing up some stuff with S
% Calculate SNR each channel across all tasks

% Save S 
if ~exist(fullfile(dataFolder, 'processed_', filename), 'dir')
    mkdir(fullfile(dataFolder, 'processed', filename));
end

%% Formatting data from each file
% Load in nev and out for given file
this_task = 'mdir1';
tic
load(sprintf('%s/nev_out/%s_%s.mat',dataFolder,filename,this_task))
toc

tic
tbl = format_dataTable(nev,out,this_task);
toc

% %%
% 
% f = figure;
% for t=1:height(tbl)
%     if tbl.result(t)=="CORRECT"
%         [th,rho] = cart2pol(tbl.eyePos{t}(1,:),tbl.eyePos{t}(2,:));
%         plot(rho,'k-','linewidth',1)
%     end
%     hold on;
%     prettyFig;
% 
%     blah = 1;
% end
% 
% %%
% if ~exist(fullfile(dataFolder, 'processed', filename), 'dir')
%     mkdir(fullfile(dataFolder, 'processed', filename));
% end
% 
% S = struct();
% 
% S.experimenter = EXPERIMENTER;
% S.monkey = MONKEY;
% S.session = SESSION;
% S.date_time = this_sess.date_time;
% S.grid_hole = this_sess.gridHole;
% S.recordDepth_mm = this_sess.recordDepth_mm;
% S.gtHeight_mm = this_sess.gtHeight_mm;
% S.probeType = this_sess.probeType;
% S.mapFile_name = this_sess.mapFile_name;
% S.firstContactdist_mm = this_sess.firstContactDist_mm;
% 
% %% PURSUIT PROCESS
% % First pass: Gather all unique column names
% allVarNames = {};
% for s = 1:this_sess.purs_num
%     this_task = sprintf('purs%d', s);
%     load(sprintf('%s/nev_dat/%s_%s.mat', dataFolder, filename, this_task))
% 
%     dat2 = dat(cellfun(@(q) ~isempty(q), {dat.nsTime}.', 'uni', 1));
%     if length(dat) ~= length(dat2)
%         fprintf('\nNEED TO REDO THIS LATER FOR NEV PAUSE\n')
%         dat = dat2;
%     end
% 
%     tbl = mayoFormat_trialDataTable(nev, dat, this_task);
% 
%     % Remove 'STIM8_ON' if it exists
%     if ismember('STIM8_ON', tbl.Properties.VariableNames)
%         tbl.STIM8_ON = [];
%     end
% 
%     % Update the list of all variable names
%     allVarNames = union(allVarNames, tbl.Properties.VariableNames, 'stable'); % 'stable' preserves order
% end
% 
% % Second pass: Standardize and concatenate tables
% purs_tbl = table(); % Initialize an empty table
% for s = 1:this_sess.purs_num
%     this_task = sprintf('purs%d', s);
%     load(sprintf('%s/nev_dat/%s_%s.mat', dataFolder, filename, this_task))
% 
%     dat2 = dat(cellfun(@(q) ~isempty(q), {dat.nsTime}.', 'uni', 1));
%     if length(dat) ~= length(dat2)
%         fprintf('\nNEED TO REDO THIS LATER FOR NEV PAUSE\n')
%         dat = dat2;
%     end
% 
%     tbl = mayoFormat_trialDataTable(nev, dat, this_task);
% 
%     % Remove 'STIM8_ON' if it exists
%     if ismember('STIM8_ON', tbl.Properties.VariableNames)
%         tbl.STIM8_ON = [];
%     end
% 
%     % Standardize the table to include all columns in `allVarNames`
%     tbl = standardizeTableColumns(tbl, allVarNames);
% 
%     % Concatenate the tables
%     purs_tbl = [purs_tbl; tbl];
% end
% 
% % List of columns to move to the end
% columnsToMove = {'eyeRaw', 'eyePos', 'eyeVel', 'eyeAcc'};
% 
% % Identify columns that exist in the final table
% columnsToMove = intersect(columnsToMove, purs_tbl.Properties.VariableNames, 'stable');
% 
% % Find the remaining columns (i.e., those not in columnsToMove)
% remainingColumns = setdiff(purs_tbl.Properties.VariableNames, columnsToMove, 'stable');
% 
% % Rearrange the columns: remaining first, then columnsToMove at the end
% purs_tbl = purs_tbl(:, [remainingColumns, columnsToMove]);
% 
% S.params = dat(1).params.block; 
% S.PURS = purs_tbl;
% 
% % Define the full file path using sprintf
% filePath = fullfile(dataFolder, 'processed', filename, sprintf('%s.mat', filename));
% 
% % Save the structure S to the specified file
% save(filePath, 'S');
% 
% %% PURSUIT PLOTS
% 
% tbl = S.PURS;
% 
% f1 = figure;
% f1.Position = [100 100 1700 600];
% tl = tiledlayout(2,3);
% tl.TileSpacing = 'tight';
% tl.Padding = 'loose';
% 
% % pie chart of trial outcomes
% results = sort(tbl.trialOutcome);
% 
% nexttile([2,1])
% pieChart_trialOutcomes(results)
% 
% % percent correct by condition
% angs = sort(unique(tbl.angle)); angs = [angs; angs(1)];
% colors = cellfun(@(q) radial_colormap(q), num2cell(angs), 'uni', 0);
% 
% jumps = sort(unique(tbl.jump));
% speeds = sort(unique(tbl.pursuitSpeed));
% perCorrects = zeros(length(angs),length(jumps),length(speeds));
% for a=1:length(angs)
%     for j=1:length(jumps)
%         for s=1:length(speeds)
%             perCorrects(a,j,s) = sum(tbl.trialOutcome=='CORRECT' & tbl.angle==angs(a) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s))/sum(tbl.angle==angs(a) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s));
%         end
%     end
% end
% 
% angs_rad = angs * (pi / 180);
% ls = {'-','--',':'};
% fa = {1, 0.75, 0.5};
% 
% for s=1:length(speeds)
%     nexttile([1,1])
%     l = gobjects(1, length(jumps)); % Preallocate array of graphics objects
%     for j=1:length(jumps)
%         this_line = perCorrects(:,j,s).*100;
% 
%         l(j) = polarplot(angs_rad,this_line, 'Color', 'black', 'LineStyle', ls{j}, 'LineWidth', 1.5);
%         hold on;
%         % Plot each point with its specific color
%         for a = 1:length(angs)
%             polarplot(angs_rad(a), this_line(a), 'o', ...
%                 'MarkerSize', 10, 'MarkerFaceColor', colors{a}, 'MarkerEdgeColor', colors{a});
%         end
%     end
%     ax = gca; % Get the polar axes
%     ax.RLim = [0 100]; % Set radial limits
%     ax.ThetaTick = 0:30:180; % Set angular ticks (degrees)
%     ax.RTick = [0 20 40 60 80 100]; % Set radial ticks
% 
%     if s==1
%         ll = legend(l, string(jumps)');  % Create the legend
%         title(ll, 'Jump Values');
%     end
% 
%     % Add labels for the axes
%     % Radial label
%     text(0, 40, '% correct', 'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', 'FontSize', 12);
% 
%     prettyFig;
%     title(sprintf('speed = %d deg/s',speeds(s)))
% end
% 
% nexttile([1,2])
% 
% % Group by fixDuration and split the trialOutcome values
% [G, fixDurValues] = findgroups(tbl.fixDuration);
% groupedTrialOutcomes = splitapply(@(x){x}, tbl.trialOutcome, G);
% 
% % Calculate the percentage of correct trials per group
% perCorrect_perGroup = cellfun(@(q) 100*(sum(q=='CORRECT')/length(q)), groupedTrialOutcomes, 'uni', 1);
% 
% % Plot the percentage of correct trials
% plot(fixDurValues, perCorrect_perGroup, 'k.-', 'LineWidth', 2, 'MarkerSize', 10)
% xlabel('Fixation duration prior to target motion onset (ms)')
% ylabel('Percent correct (%)')
% 
% prettyFig;
% % Annotate each point with the count of values in each group
% for i = 1:length(groupedTrialOutcomes)
%     % Get the count of trials in the current group
%     count = length(groupedTrialOutcomes{i});
% 
%     % Place the text annotation above the data point
%     text(fixDurValues(i), perCorrect_perGroup(i) + 4, num2str(count), ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10)
% end
% 
% exportgraphics(f1,fullfile(dataFolder, 'processed', filename, 'purs_behav_trialOutcomes.png'),'Resolution',300)
% 
% %%
% f2 = figure;
% f2.Position = [100 100 1500 900];
% tl = tiledlayout(6,3);
% tl.TileSpacing = 'compact';
% tl.Padding = 'loose';
% 
% % eye traces 
% preint = 25;
% postint = 210;
% rt = dat(1).params.block.reactionTime;
% clc
% cnt = 1;
% for s=1:length(speeds)
%     angles_eachJump = cell(1,length(jumps));
%     for j=1:length(jumps)
%         ax(cnt) = nexttile([2,1]); 
%         cnt = cnt + 1;
% 
%         eyeTraces = tbl.eyeVel(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         if ismember('PURSUIT_TARG', tbl.Properties.VariableNames)
%             targetOnsets = tbl.PURSUIT_TARG(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         else
%             targetOnsets = tbl.TARG_ON(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         end
%         csFlags = tbl.csFlag(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         angles = tbl.angle(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
% 
%         angles_eachJump{j} = angles(~logical(csFlags));
% 
%         % percent removed, including microsaccades and catch-up saccades
%         percentRemoved1 = 100*(sum(tbl.csFlag~=0 & tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s))/sum(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s)));
%         percentRemoved2 = 100*(sum(tbl.csFlag==1 & tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s))/sum(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s)));
% 
%         aligned = cellfun(@(q, w) q(:, w - preint:w + postint), eyeTraces, num2cell(targetOnsets), 'uni', 0);
% 
%         x = (1:length(aligned{1})) - preint;
%         ylim([0 max(speeds)*1.8])
% 
%         yLimits = ylim; % Get current y-axis limits
%         l(4) = fill([0 rt rt 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
% 
%         yline(speeds(s), 'r--')
%         hold on;
%         for t = 1:length(eyeTraces)
%             this_trl = aligned{t};
%             [~, rho] = cart2pol(this_trl(1, :), this_trl(2, :));
%             if csFlags(t) == 0
%                 plot(x, rho, 'k-')
%             else
%                 continue
%             end
%             hold on;
%         end
% 
%         xlim([-preint postint])
%         ylim([0 max(speeds)*1.8])
% 
%         if s == 1
%             title(sprintf('Jump = %d', jumps(j)))
%         end
%         if j==1
%             ylabel({'radial eye velocity', '(deg/s)'});
%         elseif j == 3
%             % Move ylabel to the right side of the plot
%             yl = ylabel(sprintf('Speed = %d deg/s', speeds(s)));
%             yl.Rotation = -90; % Rotate inward
%             yl.VerticalAlignment = 'top';
%             yl.HorizontalAlignment = 'center';
%             yl.Position = [max(xlim(ax(cnt - 1))) + 15, mean(ylim(ax(cnt - 1))), 0]; % Adjust position
%             yl.FontWeight = 'bold'; % Make the label bold
%         end
% 
%         % Add text box at the top-left corner of the current tile
%         text(ax(cnt-1), -preint + 10, max(speeds)*1.8-2, sprintf('%0.1f%% trials removed \n (%0.1f%% CS)', percentRemoved1, percentRemoved2), ...
%             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% 
%         xlabel('time aligned to target motion onset (ms)')
%         prettyFig;
%     end
% 
%     for jj=1:length(jumps)
%         ax(cnt) = nexttile([1,1]); 
%         cnt = cnt + 1;
%         % Define possible directions (0:45:315)
%         possibleDirections = 0:45:315;
% 
% 
%         % Count occurrences of each direction
%         counts = histcounts(angles_eachJump{jj}, [possibleDirections - 22.5, 337.5 + 22.5]);
% 
%         % Create the bar plot
%         bar(possibleDirections, counts, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
%         xticks(possibleDirections);
%         xlabel('target angle (deg)');
%         title(sprintf('N = %d', length(angles_eachJump{jj})),'FontWeight','normal')
%         if jj==1
%             ylabel('Count');
%         end
%         yline(height(tbl(tbl.trialOutcome=='CORRECT' & tbl.angle==angs(1) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s),:)),'-','Color',[0.15 0.15 0.15])
%         ylim([0 height(tbl(tbl.trialOutcome=='CORRECT' & tbl.angle==angs(1) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s),:))])
%         prettyFig;
% 
%         blah = 1;
%     end
% end
% 
% title(tl,sprintf('%s_purs',filename),'fontsize',18,'interpreter','latex')
% subtitle(tl, sprintf('(# of blocks = %d, step-ramp duration = %d)',this_sess.purs_num,S.params.reactionTime))
% 
% exportgraphics(f2,fullfile(dataFolder, 'processed', filename, 'purs_behav_csTrialDist.png'),'Resolution',300)
% 
% %%
% f3 = figure;
% f3.Position = [100 100 1500 900];
% tl = tiledlayout(6,3);
% tl.TileSpacing = 'compact';
% tl.Padding = 'loose';
% 
% % eye traces 
% preint = 50;
% postint = 210;
% rt = dat(1).params.block.reactionTime;
% 
% cnt = 1;
% for s=1:length(speeds)
%     angles_eachJump = cell(1,length(jumps));
%     for j=1:length(jumps)
%         ax(cnt) = nexttile([2,1]); 
%         cnt = cnt + 1;
% 
%         eyeTraces = tbl.eyeVel(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         if ismember('PURSUIT_TARG', tbl.Properties.VariableNames)
%             targetOnsets = tbl.PURSUIT_TARG(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         else
%             targetOnsets = tbl.TARG_ON(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         end
%         csFlags = tbl.csFlag(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
%         angles = tbl.angle(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s) & tbl.msFlag == 0 & ~isnan(tbl.csFlag));
% 
%         angles_eachJump{j} = angles(~logical(csFlags));
% 
%         % percent removed, including microsaccades and catch-up saccades
%         percentRemoved1 = 100*(sum(tbl.csFlag~=0 & tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s))/sum(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s)));
%         percentRemoved2 = 100*(sum(tbl.csFlag==1 & tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s))/sum(tbl.trialOutcome == "CORRECT" & tbl.jump == jumps(j) & tbl.pursuitSpeed == speeds(s)));
% 
%         aligned = cellfun(@(q, w) q(:, w - preint:w + postint), eyeTraces, num2cell(targetOnsets), 'uni', 0);
% 
%         x = (1:length(aligned{1})) - preint;
%         ylim([0 max(speeds)*3])
% 
%         yLimits = ylim; % Get current y-axis limits
%         l(4) = fill([0 rt rt 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
% 
%         yline(speeds(s), 'r--')
%         hold on;
%         for t = 1:length(eyeTraces)
%             this_trl = aligned{t};
%             [~, rho] = cart2pol(this_trl(1, :), this_trl(2, :));
%             if csFlags(t) == 0
%                 plot(x, rho, 'k-')
%             else
%                 plot(x, rho, 'r-')
%             end
%             hold on;
%         end
% 
%         xlim([-preint postint])
%         ylim([0 max(speeds)*3])
% 
%         if s == 1
%             title(sprintf('Jump = %d', jumps(j)))
%         end
%         if j==1
%             ylabel({'radial eye velocity', '(deg/s)'});
%         elseif j == 3
%             % Move ylabel to the right side of the plot
%             yl = ylabel(sprintf('Speed = %d deg/s', speeds(s)));
%             yl.Rotation = -90; % Rotate inward
%             yl.VerticalAlignment = 'top';
%             yl.HorizontalAlignment = 'center';
%             yl.Position = [max(xlim(ax(cnt - 1))) + 15, mean(ylim(ax(cnt - 1))), 0]; % Adjust position
%             yl.FontWeight = 'bold'; % Make the label bold
%         end
% 
%         % Add text box at the top-left corner of the current tile
%         text(ax(cnt-1), -preint + 10, max(speeds)*1.8-2, sprintf('%0.1f%% trials removed \n (%0.1f%% CS)', percentRemoved1, percentRemoved2), ...
%             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
% 
%         xlabel('time aligned to target motion onset (ms)')
%         prettyFig;
%     end
% 
%     for jj=1:length(jumps)
%         ax(cnt) = nexttile([1,1]); 
%         cnt = cnt + 1;
%         % Define possible directions (0:45:315)
%         possibleDirections = 0:45:315;
% 
% 
%         % Count occurrences of each direction
%         counts = histcounts(angles_eachJump{jj}, [possibleDirections - 22.5, 337.5 + 22.5]);
% 
%         % Create the bar plot
%         bar(possibleDirections, counts, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
%         xticks(possibleDirections);
%         xlabel('target angle (deg)');
%         title(sprintf('N = %d', length(angles_eachJump{jj})),'FontWeight','normal')
%         if jj==1
%             ylabel('Count');
%         end
%         yline(height(tbl(tbl.trialOutcome=='CORRECT' & tbl.angle==angs(1) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s),:)),'-','Color',[0.15 0.15 0.15])
%         ylim([0 height(tbl(tbl.trialOutcome=='CORRECT' & tbl.angle==angs(1) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s),:))])
%         prettyFig;
% 
%         blah = 1;
%     end
% end
% 
% title(tl,sprintf('%s_purs',filename),'fontsize',18,'interpreter','latex')
% subtitle(tl, sprintf('(# of blocks = %d, step-ramp duration = %d)',this_sess.purs_num,S.params.reactionTime))
% 
% exportgraphics(f3,fullfile(dataFolder, 'processed', filename, 'purs_behav_csTrialDist_withCS.png'),'Resolution',300)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
