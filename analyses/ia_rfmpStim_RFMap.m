function ia_rfmpStim_RFMap(data, varargin)
% Receptive field map of mean firing rate for the rfsa* RF-mapping stimuli
% (stimPos/STIM_ON), restricted to only the RF flashes that occurred in a
% given time window relative to an alignment event (e.g. the first TARG_ON),
% so flashes that happen later in the trial are excluded entirely. Produces
% the same sliding-window-over-time heatmap as ia_rfMaps, via
% format_tableToRFMap + heatMap_rfOverTime.

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
nTrialsTotal = height(thisTbl);

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

    title(tl, {
        sprintf('%s --- %s --- unit %d (channel %d)', S.sess_name, probe_label, clust, chans(u));
        sprintf('RF flashes restricted to [%d, %d] ms relative to %s', TIME_BIN(1), TIME_BIN(2), alignStr)
        }, 'fontsize',16,'interpreter','none')

    annotation('textbox', [0.77 0.89 0.2 0.1], ...
               'String', sprintf('N = %d trials\n%d flashes', height(thisTbl), nStimuli), ...
               'FontSize', 14, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'right');

    if ~isempty(FIG_PATH)
        if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end
        if SAVE_PDF
            savebigPDF(fig, fullfile(FIG_PATH, sprintf('%s_clust%04d_chan%03d.pdf', probe_label, clust, chans(u))));
        else
            savebigPNG(fig, fullfile(FIG_PATH, sprintf('%s_clust%04d_chan%03d.png', probe_label, clust, chans(u))));
        end
    end
    fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE----',PROBE_INDEX, clust))
end
fprintf('\n------------------------------\n')

end
