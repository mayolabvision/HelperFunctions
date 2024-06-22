function sa_pstheyesKKN(datafolder, datafile, alignType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for extracting and preprocessing eye traces from .ns5 files
%%%%%%%%%% INPUTS %%%%%%%%%%%
% datafolder = string/char of full path where data is stored
% datafile = datafile name
% alignType = 'stim' or 'saccade'

%%%%%%%%%% EXAMPLE %%%%%%%%%%%
% e.g. datafolder = '/Users/kendranoneman/Projects/mayo/HelperFunctions/process_recordings/example_data';
%      datafile = 'sb01pursA65650026';
% tbl = qa_dirmemKKN(datafolder,datafile,'stim')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default parameters
sigma = 10; % spike density function 
if isequal(alignType,'stim')
    preint = 200; % time before t=0 to plot
    postint = 2500; % time after t=0 to plot
else
    preint = 500; % time before t=0 to plot
    postint = 450; % time after t=0 to plot
end
num_angGroups = 8; % how many groups to bin directions in (2, 4, or 8)

% Process data in with functions in nevutils and helperfunctions
tbl = extract_mayoData(datafolder, datafile);

% Calculate pursuit onset
[eyeT,eyeR] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), tbl.eyeVel(tbl.trialOutcome=="CORRECT"), 'uni', 0);
purOnsets = cellfun(@(q,r) detect_pursuitOnset(q,r,0,300,0), eyeR, num2cell(tbl.TARG_ON(tbl.trialOutcome=="CORRECT")), 'uni', 1);
tbl.pursuitOnset(tbl.trialOutcome=="CORRECT") = purOnsets;

%% Seperating by condition
% Sort and get unique angles
angs = sort(unique(tbl.angle));

% Calculate the angular range for each group
start_angles = mod(-45 + (0:num_angGroups-1) * 360 / num_angGroups, 360);
end_angles = mod(-45 + (1:num_angGroups) * 360 / num_angGroups, 360);
ranges = [start_angles', end_angles'];

% Initialize cell array to store groups
ang_groups = cell(num_angGroups, 1);

% Iterate over each range
for i = 1:num_angGroups
    % Find indices of values within the current range
    if ranges(i, 1) < ranges(i, 2)
        indices = find(angs >= ranges(i, 1) & angs < ranges(i, 2));
    else
        indices = find(angs >= ranges(i, 1) | angs < ranges(i, 2));
    end
    % Store values within the current range
    ang_groups{i} = angs(indices);
end

%% Plot PSTH
% Create figure
fig = figure;
fig.Position = [100 100 1500 800];

% Create tiled layout based on number of angular groups
if num_angGroups == 2
    tl = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    positions = [3, 1];
    centPos = 2;
elseif num_angGroups == 4 || num_angGroups == 8
    tl = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    if num_angGroups == 4
        positions = [6, 2, 4, 8];
    else
        positions = [9, 6, 3, 2, 1, 4, 7, 8];
    end
    centPos = 5;
end

% Initialize an array to store the axes handles
axes_handles = gobjects(1, num_angGroups);

% Plot data in each subplot
for a = 1:num_angGroups
    thesetrials = tbl(tbl.trialOutcome == "CORRECT" & ismember(tbl.angle, ang_groups{a}), :);

    % Create the subplot
    ax = nexttile(tl, positions(a));
    
    % Determine spike times based on alignment type
    if isequal(alignType, 'stim')
        sptimes = cellfun(@(s, d) s - d(1), thesetrials.spks, num2cell(thesetrials.TARG_ON), 'uni', 0);
        xlab = 'Time Aligned to Stim Onset (ms)';
    elseif isequal(alignType, 'saccade')
        sptimes = cellfun(@(s, a) s - a(1), thesetrials.spks, num2cell(thesetrials.pursuitOnset), 'uni', 0);
        xlab = 'Time Aligned to Saccade Onset (ms)';
    end

    % Plot raster plot and SDF
    raster_sdf(sptimes, [-preint, postint], sigma)
    prettyFig;
    
    % Set title for subplot
    title(sprintf('[%d-%d] deg', ang_groups{a}(1), ang_groups{a}(end)), 'fontsize', 16, 'fontweight', 'bold')
    
    % Store axes handle
    axes_handles(a) = ax;
end

% Create subplot for center circle, can replace with tuning polar plot
nexttile(tl, centPos);
axx = gca;
axx.Visible = 'off';
hold on;

% Plot unit circle
t = linspace(0, 2*pi, 100);
x = cos(t);
y = sin(t);
plot(x, y, 'k');  % Plot unit circle
axis equal;
hold off;

linkaxes(axes_handles, 'y')
xlabel(tl, xlab, 'FontSize', 18)
title(tl, sprintf('%s (%s-aligned)', datafile, alignType), 'fontsize', 20, 'Interpreter', 'none')

% Save figure
saveas(fig, sprintf('%s/%s_rasterSDF_%s.png', datafolder, datafile, alignType))

%% Plot EYES
% Create figure
fig2 = figure;
fig2.Position = [100 100 1500 800];

% Create tiled layout based on number of angular groups
if num_angGroups == 2
    tl = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    positions = [3, 1];
    centPos = 2;
elseif num_angGroups == 4 || num_angGroups == 8
    tl = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    if num_angGroups == 4
        positions = [6, 2, 4, 8];
    else
        positions = [9, 6, 3, 2, 1, 4, 7, 8];
    end
    centPos = 5;
end

% Initialize an array to store the axes handles
axes_handles = gobjects(1, num_angGroups);

% Plot data in each subplot
for a = 1:num_angGroups
    thesetrials = tbl(tbl.trialOutcome == "CORRECT" & ismember(tbl.angle, ang_groups{a}), :);

    % Create the subplot
    ax = nexttile(tl, positions(a));
    
    % Determine spike times based on alignment type
    if isequal(alignType, 'stim')
        xlab = 'Time Aligned to Stim Onset (ms)';
    elseif isequal(alignType, 'saccade')
        xlab = 'Time Aligned to Saccade Onset (ms)';
    end

    
    for t = 1:height(thesetrials)
        x = (1:length(thesetrials.eyePos{t}(1,:)))-thesetrials.TARG_ON(t);
        plot(x,thesetrials.eyePos{t}(1,:),'r-','LineWidth',1)
        hold on
        plot(x,thesetrials.eyePos{t}(2,:),'b-','LineWidth',1)
    end   

    xline(0,'k--','linewidth',1)
    xlim([-preint,postint])
    ylim([-25,25])
    prettyFig;

    %[mn,~,yu,yl] = sem_errorbar(eyevels);
    %fill([x fliplr(x)], [yu fliplr(yl)], sha, 'linestyle','none','FaceAlpha',0.5)

    % Set title for subplot
    title(sprintf('[%d-%d] deg', ang_groups{a}(1), ang_groups{a}(end)), 'fontsize', 16, 'fontweight', 'bold')
    
    % Store axes handle
    axes_handles(a) = ax;
end

% Create subplot for center circle, can replace with tuning polar plot
nexttile(tl, centPos);
axx = gca;
axx.Visible = 'off';
hold on;

% Plot unit circle
t = linspace(0, 2*pi, 100);
x = cos(t);
y = sin(t);
plot(x, y, 'k');  % Plot unit circle
axis equal;
hold off;

linkaxes(axes_handles, 'y')
xlabel(tl, xlab, 'FontSize', 18)
title(tl, sprintf('%s (%s-aligned)', datafile, alignType), 'fontsize', 20, 'Interpreter', 'none')

% Save figure
saveas(fig2, sprintf('%s/%s_eyePos_%s.png', datafolder, datafile, alignType))

end
