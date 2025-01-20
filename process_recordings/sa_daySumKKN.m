clear
clc

addpath(genpath('/Users/kendranoneman/Packages')) % add nevUtils and HelperFunctions to path
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions')) 

DATAFOLDER = '/Users/kendranoneman/OneDrive/DATA';
EXPERIMENTER = 'kendra';
MONKEY = 'scrappy';
SESSION = '0076a';

recordings = readtable(sprintf('%s/RECORDING_INFO.csv',DATAFOLDER));
recordings.experimenter = categorical(recordings.experimenter);
recordings.monkey = categorical(recordings.monkey);
recordings.session_depth = categorical(recordings.session_depth);
recordings.mapFile_name = categorical(recordings.mapFile_name );
recordings.gridHole = cellfun(@(q) eval(strrep(strrep(q, '“', '"'), '”', '"')), recordings.gridHole, 'uni', 0);

this_sess = recordings(recordings.experimenter==EXPERIMENTER & recordings.monkey==MONKEY & recordings.session_depth==SESSION,:);
filename = sprintf('%s_%s_%s',EXPERIMENTER,MONKEY,SESSION);

%%

taskTypes = {'purs'};
numTasks = [this_sess.purs_num];

tasks = arrayfun(@(x) sprintf('%s%d', taskTypes{1}, x), 1:numTasks(1), 'UniformOutput', false);

for f=1:length(tasks)
    this_task = tasks{f}

    if exist(sprintf('%s/nev_dat/%s_%s.mat',DATAFOLDER,filename,this_task), 'file')==2
        fprintf('nev_dat for %s already exists.\n', this_task)
    else
        fprintf('generating nev_dat for %s.\n', this_task)
        nevname = sprintf('%s/raw/%s_%s.nev',DATAFOLDER,filename,this_task);

        nev = readNEV(nevname);
        dat = nev2dat(nevname,'readNS5',true,'convertEyes',true,'allowNevPause',true,'include_0_255',true);

        save(sprintf('%s/nev_dat/%s_%s.mat',DATAFOLDER,filename,this_task), 'nev','dat')
    end
end

%%
this_task = 'purs1';
load(sprintf('%s/nev_dat/%s_%s.mat',DATAFOLDER,filename,this_task))

%% PURSUIT PROCESS
this_task = 'purs1';
tbl = mayoFormat_trialDataTable(nev,dat,this_task);


%% PURSUIT PLOTS
f1 = figure;
f1.Position = [100 100 1400 900];
t = tiledlayout(2,3);
t.TileSpacing = 'loose';
t.Padding = 'compact';

nexttile
pieChart_trialOutcomes(tbl.trialOutcome)


% percent correct by condition
angs = sort(unique(tbl.angle)); angs = [angs; angs(1)];
colors = cellfun(@(q) radial_colormap(q), num2cell(angs), 'uni', 0);

jumps = sort(unique(tbl.jump));
speeds = sort(unique(tbl.pursuitSpeed));
perCorrects = zeros(length(angs),length(jumps),length(speeds));
for a=1:length(angs)
    for j=1:length(jumps)
        for s=1:length(speeds)
            perCorrects(a,j,s) = sum(tbl.trialOutcome=='CORRECT' & tbl.angle==angs(a) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s))/sum(tbl.angle==angs(a) & tbl.jump==jumps(j) & tbl.pursuitSpeed==speeds(s));
        end
    end
end

angs_rad = angs * (pi / 180);

for s=1:length(speeds)
    nexttile
    for j=1:length(jumps)
        this_line = perCorrects(:,j,s).*100;

        polarplot(angs_rad,this_line, 'k-', 'LineWidth', 1.5);
        hold on;
        % Plot each point with its specific color
        for a = 1:length(angs)
            polarplot(angs_rad(a), this_line(a), 'o', ...
                'MarkerSize', 10, 'MarkerFaceColor', colors{a}, 'MarkerEdgeColor', colors{a});
        end
    end
    ax = gca; % Get the polar axes
    ax.RLim = [0 100]; % Set radial limits
    ax.ThetaTick = 0:30:180; % Set angular ticks (degrees)
    ax.RTick = [0 20 40 60 80 100]; % Set radial ticks
    
    % Add labels for the axes
    % Radial label
    text(0, 40, '% Correct', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12);
    
    prettyFig;
    title(sprintf('speed = %d deg/s',speeds(s)))
end

%%
% Customize axes
ax = gca; % Get the polar axes
ax.RLim = [0 100]; % Set radial limits
ax.ThetaTick = 0:30:180; % Set angular ticks (degrees)
ax.RTick = [0 20 40 60 80 100]; % Set radial ticks

% Add labels for the axes
% Radial label
text(0, 40, '% Correct', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 12);

prettyFig;


% distribution of RT by jumpSize






















