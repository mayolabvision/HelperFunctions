function results = ia_rfmpStim_RFMapFit(data, varargin)
% Same rfsa* RF map as ia_rfmpStim_RFMap, but per cluster also:
%   (1) tests whether the cluster has a statistically significant receptive
%       field at the peak-response time bin, via a nested-model F-test
%       comparing the fitted Gaussian (6 params) against a flat/no-RF null
%       (1 param: the grand mean), over the grid's cell means
%   (2) finds the RF latency -- the 50ms sliding-window time bin (relative
%       to each RF flash's own onset) with the largest deviation from that
%       bin's own grid-wide mean. Excitation and inhibition are allowed to
%       peak at different latencies, so each direction gets its own
%       independently-detected bin.
%   (3) fits a 2D Gaussian to the mean firing-rate map at each direction's
%       own peak bin to get the RF's peak position and extent
%       (fit_gaussianRF.m) -- a cluster can have an excitatory RF (region
%       of significantly elevated firing), an inhibitory RF (region of
%       significantly suppressed firing), both (each at its own latency),
%       or neither
%
% A candidate RF only counts if (a) the F-test at that direction's peak bin
% is significant, (b) the fit's R^2 clears MIN_RSQ, (c) its amplitude is
% physiologically plausible (|amp| <= MAX_AMP -- catches outlier-driven
% spurious fits), and (d) the fit didn't park at an optimizer bound instead
% of converging (fit_gaussianRF.m's gof.hitBound -- amp ~ 0, or sigma pinned
% to its generous upper cap). A peak bin's tile is only highlighted in the
% plot when a real RF is found there.
%
% (An earlier version used a one-way ANOVA of raw per-flash firing rate
% across all 144 stimulus positions to test significance, and rejected any
% fit whose center (x0/y0) landed outside the tested stimulus range. Both
% were found to reject real RFs: the ANOVA is underpowered for smooth,
% modest-firing-rate bumps once diluted across 143 degrees of freedom of
% sparse single-flash spike counts, and a real RF's center can legitimately
% sit just past the edge of the tested grid, an eccentric RF the display
% couldn't fully capture. The F-test targets the actual question -- does a
% smooth 2D bump explain the data better than no RF at all -- with far fewer
% effective parameters, and MAX_AMP targets implausible-amplitude artifacts
% directly instead of inferring them indirectly from fit location.)
%
% All of this uses only the RF flashes restricted to TIME_BIN relative to
% ALIGN_TO (default [-400,0] ms relative to the first TARG_ON), i.e. before
% anything task-related has happened. Produces the same sliding-window
% heatmap as ia_rfmpStim_RFMap, with each RF's peak bin bordered red
% (excitatory) or blue (inhibitory) -- black if both land on the same bin
% -- and the fitted Gaussian(s)' half-max contour drawn on top in matching colors.

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
addParameter(p, 'ALPHA', 0.05, @isnumeric); % significance threshold for the RF F-test
addParameter(p, 'MIN_RSQ', 0.2, @isnumeric); % minimum Gaussian fit R^2 for a candidate RF to count
addParameter(p, 'MAX_AMP', 500, @isnumeric); % max plausible |amplitude| (Hz) for a candidate RF to count
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
ALPHA = p.Results.ALPHA;
MIN_RSQ = p.Results.MIN_RSQ;
MAX_AMP = p.Results.MAX_AMP;
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

[thisTbl, alignStr, nTrialsTotal, nStimuli] = prep_rfsaTable(S, ALIGN_TO, ALIGN_IND, TIME_BIN, ANGLE);

prb_name = sprintf('spiketimes_%d', PROBE_INDEX);
if ~ismember(prb_name, thisTbl.Properties.VariableNames)
    error('Column ''%s'' not found -- these rfsa tables don''t have spike times merged in for PROBE_INDEX %d.', prb_name, PROBE_INDEX);
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
% ia_rfMaps' FIG_PATH2 (<hardware_config>_<probe_label>/rfsa_RFfits/<bin tag>),
% so re-running this analysis with different bins doesn't overwrite past runs
if ~isempty(FIG_PATH)
    binTag = sprintf('%s_%dto%dms', ALIGN_TO, TIME_BIN(1), TIME_BIN(2));
    FIG_PATH = fullfile(FIG_PATH, sprintf('%s_%s', hardware_config, probe_label), 'rfsa_RFfits', binTag);
end

results = table();

for u = 1:numel(clusts)
    clust = clusts(u);

    [frsRaw, bin_edges, xvals, yvals] = format_tableToRFMap(thisTbl, 'PROBE_INDEX', PROBE_INDEX, ...
        'UNITS', (clust+1), 'FIRST_BIN', FIRST_BIN, 'BIN_WIDTH', BIN_WIDTH, 'BIN_STEP', BIN_STEP, 'N_BINS', N_BINS);
    frsRaw = frsRaw{1}; % Ny x Nx x N_BINS cell, each cell a vector of per-flash Hz values
    meanFRs = cellfun(@mean, frsRaw); % Ny x Nx x N_BINS
    % (2) Latency: excitation and inhibition can peak at different times, so
    % find each direction's own strongest bin (largest deviation from that
    % bin's own grid-wide mean, in that direction) independently
    gridMeanPerBin = squeeze(mean(mean(meanFRs,1),2));
    devFromMean = meanFRs - reshape(gridMeanPerBin,1,1,[]);
    [~, excBinIdx] = max(squeeze(max(max(devFromMean,[],1),[],2)));
    [~, inhBinIdx] = min(squeeze(min(min(devFromMean,[],1),[],2)));

    % (1) & (3): test significance and fit a Gaussian for each direction at
    % its own peak bin (see analyzeDirection below)
    [hasExcRF, excParams, excRsq, excP] = analyzeRFDirection(meanFRs, xvals, yvals, excBinIdx, 'exc', ALPHA, MIN_RSQ, MAX_AMP);
    [hasInhRF, inhParams, inhRsq, inhP] = analyzeRFDirection(meanFRs, xvals, yvals, inhBinIdx, 'inh', ALPHA, MIN_RSQ, MAX_AMP);

    results = [results; table(clust, chans(u), ...
        excBinIdx, bin_edges{excBinIdx}(1), bin_edges{excBinIdx}(2), excP, hasExcRF, excParams(1), excParams(2), excParams(3), excParams(4), excParams(5), excRsq, ...
        inhBinIdx, bin_edges{inhBinIdx}(1), bin_edges{inhBinIdx}(2), inhP, hasInhRF, inhParams(1), inhParams(2), inhParams(3), inhParams(4), inhParams(5), inhRsq, ...
        'VariableNames', {'clust','chan', ...
        'excBinIdx','excLatencyStart','excLatencyEnd','excPValue','hasExcRF','excAmp','excX0','excY0','excSigX','excSigY','excRsq', ...
        'inhBinIdx','inhLatencyStart','inhLatencyEnd','inhPValue','hasInhRF','inhAmp','inhX0','inhY0','inhSigX','inhSigY','inhRsq'})]; %#ok<AGROW>

    if ~isempty(FIG_PATH)
        fig = figure('Visible','off');
    else
        fig = figure('Visible','on');
    end
    fig.Position = [100 100 1800 900];
    [tl, ax] = heatMap_rfOverTime(meanFRs, 'BIN_EDGES', bin_edges, 'INTERP', false, 'X_VALS', xvals, 'Y_VALS', yvals);

    % Only highlight a bin's tile if there's a real RF there. If exc and inh
    % land on the same bin, that tile gets both overlays and a black border;
    % otherwise each gets its own tile, bordered/overlaid in its own color.
    % The drawn contour is the fitted Gaussian's half-max (50%-of-peak)
    % boundary -- the standard convention for reporting RF size/extent --
    % rather than the 1 SD contour (radius = sigma * sqrt(2*ln(2))).
    theta = linspace(0, 2*pi, 100);
    hwhmFactor = sqrt(2*log(2));
    sameBin = hasExcRF && hasInhRF && excBinIdx == inhBinIdx;

    if hasExcRF
        if sameBin
            excBorderColor = 'k'; excLabel = 'EXC+INH';
        else
            excBorderColor = 'r'; excLabel = 'EXC';
        end
        set(ax(excBinIdx), 'LineWidth', 3, 'XColor', excBorderColor, 'YColor', excBorderColor, 'Box', 'on');
        title(ax(excBinIdx), sprintf('[%d , %d) ms PEAK (%s)', bin_edges{excBinIdx}(1), bin_edges{excBinIdx}(2), excLabel), ...
            'FontWeight', 'bold', 'Color', excBorderColor);
        hold(ax(excBinIdx), 'on');
        plot(ax(excBinIdx), excParams(2) + hwhmFactor*excParams(4)*cos(theta), excParams(3) + hwhmFactor*excParams(5)*sin(theta), 'r-', 'LineWidth', 2);
        plot(ax(excBinIdx), excParams(2), excParams(3), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        hold(ax(excBinIdx), 'off');
    end

    if hasInhRF
        if ~sameBin
            set(ax(inhBinIdx), 'LineWidth', 3, 'XColor', 'b', 'YColor', 'b', 'Box', 'on');
            title(ax(inhBinIdx), sprintf('[%d , %d) ms PEAK (INH)', bin_edges{inhBinIdx}(1), bin_edges{inhBinIdx}(2)), ...
                'FontWeight', 'bold', 'Color', 'b');
        end
        hold(ax(inhBinIdx), 'on');
        plot(ax(inhBinIdx), inhParams(2) + hwhmFactor*inhParams(4)*cos(theta), inhParams(3) + hwhmFactor*inhParams(5)*sin(theta), 'b-', 'LineWidth', 2);
        plot(ax(inhBinIdx), inhParams(2), inhParams(3), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
        hold(ax(inhBinIdx), 'off');
    end

    if isnan(ANGLE)
        angleTitleStr = 'all angles';
    else
        angleTitleStr = sprintf('angle = %g', ANGLE);
    end
    title(tl, {
        sprintf('%s --- %s --- unit %d (channel %d)', S.sess_name, probe_label, clust, chans(u));
        sprintf('RF flashes restricted to [%d, %d] ms relative to %s, %s', TIME_BIN(1), TIME_BIN(2), alignStr, angleTitleStr)
        }, 'fontsize',16,'interpreter','none')

    rfLines = {};
    if hasExcRF
        rfLines{end+1} = sprintf('Excitatory RF: latency = [%d, %d) ms, center = (%.0f, %.0f), half-max = (%.0f, %.0f)', ...
            bin_edges{excBinIdx}(1), bin_edges{excBinIdx}(2), excParams(2), excParams(3), hwhmFactor*excParams(4), hwhmFactor*excParams(5));
    end
    if hasInhRF
        rfLines{end+1} = sprintf('Inhibitory RF: latency = [%d, %d) ms, center = (%.0f, %.0f), half-max = (%.0f, %.0f)', ...
            bin_edges{inhBinIdx}(1), bin_edges{inhBinIdx}(2), inhParams(2), inhParams(3), hwhmFactor*inhParams(4), hwhmFactor*inhParams(5));
    end
    if isempty(rfLines)
        rfLines{1} = sprintf('no RF found (exc p = %.3g, inh p = %.3g)', excP, inhP);
    end
    annotation('textbox', [0.75 0.83 0.24 0.16], ...
               'String', sprintf('N = %d trials, %d flashes\n%s', height(thisTbl), nStimuli, strjoin(rfLines, '\n')), ...
               'FontSize', 12, ...
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
    fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE (exc: p=%.3g @[%d,%d)ms=%d, inh: p=%.3g @[%d,%d)ms=%d)----', ...
        PROBE_INDEX, clust, excP, bin_edges{excBinIdx}(1), bin_edges{excBinIdx}(2), hasExcRF, ...
        inhP, bin_edges{inhBinIdx}(1), bin_edges{inhBinIdx}(2), hasInhRF))
end
fprintf('\n------------------------------\n')

end
