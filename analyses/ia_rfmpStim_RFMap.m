function ia_rfmpStim_RFMap(data, varargin)
% ia_rfmpStim_RFMap — Sliding-window receptive field map of mean firing rate
% for the rfsa* RF-mapping stimuli (stimPos/STIM_ON), restricted to only the
% RF flashes that occurred in a given time window relative to an alignment
% event (e.g. the first TARG_ON), so flashes that happen later in the trial
% are excluded entirely. Produces the same sliding-window-over-time heatmap
% as ia_rfMaps, via format_tableToRFMap + heatMap_rfOverTime.
%
%%% Required Input: %%%
%   data  -  Either:
%              - a char/string path to a session .mat file containing a
%                struct 'S' (the file is loaded internally), or
%              - the already-loaded session struct S itself.
%            S must contain:
%              - one or more fields whose name contains 'rfsa' (e.g. rfsa1,
%                rfsa2, ...), each with a '.tbl' table of trials. Tables
%                across multiple rfsa fields are combined automatically.
%                Each table needs (at minimum): result, angle, stimPos,
%                STIM_ON, an alignment-event column (default TARG_ON), and
%                params (used to convert stimPos from pixels to degrees).
%              - S.sorting: a struct array, one entry per probe, with
%                .probe_index and .clusters (a table with cluster_id, and
%                optionally best_channel/probe_label/hardware_config).
%              - a spiketimes_<PROBE_INDEX> column in the rfsa table(s),
%                merged in during spike-sorting (this function does not sort
%                spikes itself -- if that column is missing you'll get a
%                clear error telling you so).
%
%%% Example Usage: %%%
%   1. Load S into the workspace
%   load('/Volumes/lab_NHPdata-processed/kendra_scrappy_0177a_g0/tables/kendra_scrappy_0177a_g0-1b86ffd01dc9c878.mat');
%
%   2. Run function to create RF map
%   ia_rfmpStim_RFMap(S, 'PROBE_INDEX', 1, 'CLUSTER', 0);
%
%   % or pass the file path directly instead of loading it yourself first:
%   ia_rfmpStim_RFMap('PATH/TO/DATA/kendra_scrappy_0177a_g0-1b86ffd01dc9c878.mat', ...
%       'PROBE_INDEX', 1, 'CLUSTER', 0);
%
%%% Optional Name-Value Pair Arguments: %%%
%   'PROBE_INDEX' -  Which probe to use (default 1). This is MATLAB's
%                    1-indexed probe number as stored in S.sorting.probe_index,
%                    so it can never be lower than 1. Check
%                    unique([S.sorting.probe_index]) if you're not sure how
%                    many probes a session has or what indices they use.
%   'CLUSTER'     -  Which cluster ID to plot (default: [], i.e. unset).
%                    Cluster IDs are 0-indexed (as output by the spike
%                    sorter), so CLUSTER=0 is a real, valid cluster -- not
%                    the same thing as leaving CLUSTER unspecified. If you
%                    leave CLUSTER empty, the function loops through and
%                    plots EVERY cluster recorded on PROBE_INDEX, which can
%                    be hundreds to thousands of clusters and take a long
%                    time. Specify a single CLUSTER while exploring; only
%                    leave it empty once you actually want the full batch
%                    (see JOB_ID/N_CHUNKS below for splitting that batch
%                    across a SLURM array job instead of running it all in
%                    one call).
%   'ALIGN_TO'    -  Column in the trial table to align STIM_ON to (default
%                    'TARG_ON'). Can be a column holding a single numeric
%                    time per trial (e.g. 'SACCADE'), or a column holding a
%                    cell array of multiple times per trial (e.g. 'TARG_ON',
%                    which can fire more than once) -- see ALIGN_IND below
%                    for the latter case.
%   'ALIGN_IND'   -  If ALIGN_TO is a cell-per-trial column, which value
%                    within that trial's cell to align to (default 1, i.e.
%                    the first occurrence).
%   'TIME_BIN'    -  [start end] in ms, relative to ALIGN_TO(ALIGN_IND)
%                    (default [-400 0]). Only RF flashes whose STIM_ON falls
%                    in this window are used to build the map -- e.g. the
%                    default only uses flashes from 400ms before up to the
%                    moment of the first TARG_ON, i.e. before anything
%                    task-related has happened yet.
%   'FIRST_BIN'   -  Sliding-window PSTH parameters (ms), passed straight
%   'BIN_WIDTH'      through to format_tableToRFMap: FIRST_BIN is the first
%   'BIN_STEP'       bin's start time relative to each individual flash's own
%   'N_BINS'         onset (default 0), BIN_WIDTH is each bin's width
%                    (default 50), BIN_STEP is the step between consecutive
%                    bins (default 10), and N_BINS is how many bins to
%                    compute (default 24) -- so by default this produces 24
%                    overlapping 50ms-wide bins spanning 0 to 280ms
%                    post-flash.
%   'ANGLE'       -  If specified, only use trials where thisTbl.angle
%                    equals this value (default NaN, meaning use all
%                    angles).
%   'FIG_PATH'    -  If given, saves each cluster's figure as a .png (or
%                    .pdf, see SAVE_PDF) instead of displaying it on screen,
%                    under
%                    FIG_PATH/<hardware_config>_<probe_label>/rfsa_heatmaps/<ALIGN_TO>_<TIME_BIN(1)>to<TIME_BIN(2)>ms/.
%                    Default [] (no path), which instead pops up a visible
%                    figure per cluster -- fine when CLUSTER is a single ID,
%                    but avoid leaving both FIG_PATH and CLUSTER empty at
%                    the same time, since that would try to pop up one
%                    figure window per cluster in the whole batch.
%   'SAVE_PDF'    -  If true, saves as .pdf instead of .png (default false).
%                    Only relevant when FIG_PATH is given.
%   'JOB_ID'      -  For splitting the full cluster list across a SLURM array
%   'N_CHUNKS'       job: JOB_ID is the 0-indexed task ID (e.g.
%                    $SLURM_ARRAY_TASK_ID), N_CHUNKS is the total number of
%                    tasks in the array (e.g. $SLURM_ARRAY_TASK_COUNT). Only
%                    used when CLUSTER is left empty; ignored otherwise.
%                    Default NaN for both, meaning "run every cluster on
%                    PROBE_INDEX in this one call."

p = inputParser;
addRequired(p, 'data', @(x) ischar(x) || isstruct(x));
addParameter(p, 'PROBE_INDEX', 1, @isnumeric);
addParameter(p, 'CLUSTER', [], @isnumeric);
addParameter(p, 'ALIGN_TO', 'TARG_ON', @ischar);
addParameter(p, 'ALIGN_IND', 1, @isnumeric);
addParameter(p, 'TIME_BIN', [-400 0], @isnumeric); % ms rel. to ALIGN_TO -- only STIM_ON flashes in this window are used
addParameter(p, 'FIRST_BIN', 0, @isnumeric); % sliding-window PSTH params, passed through to format_tableToRFMap
addParameter(p, 'BIN_WIDTH', 50, @isnumeric);
addParameter(p, 'BIN_STEP', 10, @isnumeric);
addParameter(p, 'N_BINS', 24, @isnumeric);
addParameter(p, 'FIG_PATH', [], @ischar);
addParameter(p, 'SAVE_PDF', false, @islogical);
addParameter(p, 'JOB_ID', NaN, @isnumeric);
addParameter(p, 'N_CHUNKS', NaN, @isnumeric);
addParameter(p, 'ANGLE', NaN, @isnumeric); % if specified, only use trials where thisTbl.angle == ANGLE

parse(p, data, varargin{:});
data = p.Results.data;
PROBE_INDEX = p.Results.PROBE_INDEX;
CLUSTER = p.Results.CLUSTER;
ALIGN_TO = p.Results.ALIGN_TO;
ALIGN_IND = p.Results.ALIGN_IND;
TIME_BIN = p.Results.TIME_BIN;
FIRST_BIN = p.Results.FIRST_BIN;
BIN_WIDTH = p.Results.BIN_WIDTH;
BIN_STEP = p.Results.BIN_STEP;
N_BINS = p.Results.N_BINS;
FIG_PATH = p.Results.FIG_PATH;
SAVE_PDF = p.Results.SAVE_PDF;
JOB_ID = p.Results.JOB_ID;
N_CHUNKS = p.Results.N_CHUNKS;
ANGLE = p.Results.ANGLE;

fprintf('\n------------------------------\n')
if ischar(data)
    [~, filename, ~] = fileparts(data);
    load(data,'S');
    fprintf(sprintf('\n----Data loaded for %s----\n',filename))
else
    S = data;
end

% Combine every task struct whose field name contains 'rfsa' (rfsa1, rfsa2, ...)
sFields = fieldnames(S);
rfsaFields = sFields(contains(sFields, 'rfsa', 'IgnoreCase', true));
if isempty(rfsaFields)
    error('No fields containing ''rfsa'' were found in the data struct.');
end
fprintf('Combining %d rfsa table(s): %s\n', numel(rfsaFields), strjoin(rfsaFields, ', '));

rfsaTbls = cell(numel(rfsaFields),1);
for i = 1:numel(rfsaFields)
    rfsaTbls{i} = S.(rfsaFields{i}).tbl;
end
thisTbl = vertcat(rfsaTbls{:});
thisTbl = thisTbl(thisTbl.result=='CORRECT',:);

if ~isnan(ANGLE)
    if ~ismember(ANGLE, thisTbl.angle)
        error('ANGLE %g not found among thisTbl.angle values: %s', ANGLE, mat2str(unique(thisTbl.angle)'));
    end
    thisTbl = thisTbl(thisTbl.angle==ANGLE,:);
    fprintf('Restricting to angle == %g (%d trials)\n', ANGLE, height(thisTbl));
end
nTrialsTotal = height(thisTbl);

% stimPos is stored in screen pixels -- convert to degrees of visual angle
screenDist = thisTbl.params(1,1).block.screenDistance;
pixPerCM = thisTbl.params(1,1).block.pixPerCM;
thisTbl.stimPos = cellfun(@(pos) pix2deg(pos, screenDist, pixPerCM), thisTbl.stimPos, 'UniformOutput', false);

prb_name = sprintf('spiketimes_%d', PROBE_INDEX);
if ~ismember(prb_name, thisTbl.Properties.VariableNames)
    error('Column ''%s'' not found -- these rfsa tables don''t have spike times merged in for PROBE_INDEX %d.', prb_name, PROBE_INDEX);
end

% ALIGN_TO must be on the same ms clock as STIM_ON. It can either hold a
% single scalar per trial (e.g. SACCADE), or a cell array of multiple
% values per trial (e.g. TARG_ON) -- in the latter case ALIGN_IND picks
% which value within that trial's cell to align to
if ~ismember(ALIGN_TO, thisTbl.Properties.VariableNames)
    error('ALIGN_TO ''%s'' is not a column of thisTbl.', ALIGN_TO);
end
alignColAll = thisTbl.(ALIGN_TO);
if isnumeric(alignColAll)
    isCellAlign = false;
elseif iscell(alignColAll)
    isCellAlign = true;
else
    error('ALIGN_TO ''%s'' must be numeric-per-trial or a cell array of numeric vectors, not %s.', ALIGN_TO, class(alignColAll));
end

% Keep only the RF flashes (stimPos/STIM_ON) that fall within TIME_BIN
% relative to ALIGN_TO(ALIGN_IND) -- flashes elsewhere in the trial (e.g.
% after the target appears) are dropped before computing firing rates.
% format_tableToRFMap expects a 'conditions' column: per trial, a 1xN cell
% array where each cell is a 1x2 [x y] stimulus position (paired by index
% with that trial's STIM_ON entries).
thisTbl.conditions = cell(height(thisTbl),1);
for t = 1:height(thisTbl)
    pos = thisTbl.stimPos{t};     % Nx2 [x y]
    onTimes = thisTbl.STIM_ON{t}; % Nx1, ms relative to trial onset

    if isCellAlign
        alignVals = thisTbl.(ALIGN_TO){t};
        if numel(alignVals) >= ALIGN_IND
            eventTime = alignVals(ALIGN_IND);
        else
            eventTime = NaN; % trial has fewer ALIGN_TO events than ALIGN_IND -- contributes nothing
        end
    else
        eventTime = thisTbl.(ALIGN_TO)(t);
    end

    relOnTimes = onTimes - eventTime;
    inBin = relOnTimes >= TIME_BIN(1) & relOnTimes < TIME_BIN(2);

    thisTbl.STIM_ON{t} = onTimes(inBin);
    thisTbl.conditions{t} = num2cell(pos(inBin,:), 2);
end

% Drop trials left with no qualifying flashes
thisTbl = thisTbl(cellfun(@numel, thisTbl.STIM_ON) > 0, :);
nStimuli = sum(cellfun(@numel, thisTbl.STIM_ON));
if isCellAlign
    alignStr = sprintf('%s(%d)', ALIGN_TO, ALIGN_IND);
else
    alignStr = ALIGN_TO;
end
fprintf('%d of %d correct trials have >=1 RF flash in [%d, %d] ms relative to %s (%d flashes total)\n', ...
    height(thisTbl), nTrialsTotal, TIME_BIN(1), TIME_BIN(2), alignStr, nStimuli);
if isempty(thisTbl)
    error('No RF flashes fall in [%d, %d] ms relative to %s -- nothing to plot.', TIME_BIN(1), TIME_BIN(2), alignStr);
end

% Identify clusters/units for this probe, same convention as ia_rfMaps
if ~isfield(S, 'sorting')
    error('S.sorting not found -- spike sorting data is required to compute firing rates.');
end
sortEntry = S.sorting([S.sorting.probe_index] == PROBE_INDEX);
if isempty(sortEntry)
    error('No S.sorting entry found for PROBE_INDEX %d.', PROBE_INDEX);
end
clusts_all = sortEntry.clusters.cluster_id;
if ismember('best_channel', sortEntry.clusters.Properties.VariableNames)
    chans = sortEntry.clusters.best_channel;
else
    chans = nan(numel(clusts_all),1);
end

if isempty(CLUSTER)
    if ~isnan(JOB_ID)
        all_units = clusts_all + 1;
        % Split into N_CHUNKS chunks as a cell array, one chunk per SLURM array task
        chunks = arrayfun(@(i) all_units(...
            floor((i-1)*numel(all_units)/N_CHUNKS)+1 : ...
            floor(i*numel(all_units)/N_CHUNKS)), ...
            1:N_CHUNKS, 'UniformOutput', false);
        ids = chunks{(JOB_ID+1)};

        clusts = clusts_all(ids);
        chans = chans(ids);
    else
        clusts = clusts_all;
    end
else
    clusts = CLUSTER;
    chans = chans(clusts_all==CLUSTER);
end
if isempty(chans)
    error('CLUSTER %d not found among recorded clusters for PROBE_INDEX %d.', CLUSTER, PROBE_INDEX);
end

probe_label = string(sortEntry.clusters.probe_label(1));
hardware_config = string(sortEntry.clusters.hardware_config(1));

% Save each ALIGN_TO/TIME_BIN combo to its own subfolder, same layout as
% ia_rfMaps' FIG_PATH2 (<hardware_config>_<probe_label>/rfsa_heatmaps/<bin tag>),
% so re-running this analysis with different bins doesn't overwrite past runs
if ~isempty(FIG_PATH)
    binTag = sprintf('%s_%dto%dms', ALIGN_TO, TIME_BIN(1), TIME_BIN(2));
    FIG_PATH = fullfile(FIG_PATH, sprintf('%s_%s', hardware_config, probe_label), 'rfsa_heatmaps', binTag);
end

for u = 1:numel(clusts)
    clust = clusts(u);

    [frs, bin_edges, xvals, yvals] = format_tableToRFMap(thisTbl, 'PROBE_INDEX', PROBE_INDEX, ...
        'UNITS', (clust+1), 'FIRST_BIN', FIRST_BIN, 'BIN_WIDTH', BIN_WIDTH, 'BIN_STEP', BIN_STEP, 'N_BINS', N_BINS);

    if ~isempty(FIG_PATH)
        fig = figure('Visible','off');
    else
        fig = figure('Visible','on');
    end
    fig.Position = [100 100 1800 900];
    tl = heatMap_rfOverTime(frs{1}, 'BIN_EDGES', bin_edges, 'INTERP', false, 'X_VALS', xvals, 'Y_VALS', yvals);

    if isnan(ANGLE)
        angleTitleStr = 'all angles';
    else
        angleTitleStr = sprintf('angle = %g', ANGLE);
    end
    title(tl, {
        sprintf('%s --- %s --- unit %d (channel %d)', S.sess_name, probe_label, clust, chans(u));
        sprintf('RF flashes restricted to [%d, %d] ms relative to %s, %s', TIME_BIN(1), TIME_BIN(2), alignStr, angleTitleStr)
        }, 'fontsize',16,'interpreter','none')

    annotation('textbox', [0.77 0.89 0.2 0.1], ...
               'String', sprintf('N = %d trials\n%d flashes', height(thisTbl), nStimuli), ...
               'FontSize', 14, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'right');

    if ~isempty(FIG_PATH)
        if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end
        fileBase = sprintf('%s_clust%04d_chan%03d', probe_label, clust, chans(u));
        if ~isnan(ANGLE)
            fileBase = sprintf('%s_ang%03d', fileBase, ANGLE);
        end
        if SAVE_PDF
            savebigPDF(fig, fullfile(FIG_PATH, [fileBase '.pdf']));
        else
            savebigPNG(fig, fullfile(FIG_PATH, [fileBase '.png']));
        end
    end
    fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE----',PROBE_INDEX, clust))
end
fprintf('\n------------------------------\n')

end
