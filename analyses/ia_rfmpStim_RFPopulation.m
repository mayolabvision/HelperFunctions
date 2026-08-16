function popResults = ia_rfmpStim_RFPopulation(data, varargin)
% For every cluster on every probe, runs the same latency-detection,
% significance-testing, and Gaussian-fit analysis as ia_rfmpStim_RFMapFit
% (reusing prep_rfsaTable.m and analyzeRFDirection.m), but does NOT generate
% or save any per-cluster RF map figures -- only a single population summary
% figure is produced, a 2 x nProbes tiled layout:
%   top row    : for each probe, every detected RF's half-max contour
%                overlaid on one plot (red = excitatory, blue = inhibitory),
%                to show the population's spatial RF coverage
%   bottom row : for each probe, a histogram of RF latencies (excitatory in
%                red, inhibitory in blue)
%
% Since the shared trial-table prep (combining rfsa tables, pix2deg
% conversion, TIME_BIN/ALIGN_TO flash restriction) is identical across every
% cluster, it's done once up front rather than once per cluster.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'data', @(x) ischar(x) || isstruct(x));
addParameter(p, 'PROBES', [], @isnumeric); % default: every probe_index in S.sorting
addParameter(p, 'ALIGN_TO', 'TARG_ON', @ischar);
addParameter(p, 'ALIGN_IND', 1, @isnumeric);
addParameter(p, 'TIME_BIN', [-400 0], @isnumeric);
addParameter(p, 'FIRST_BIN', 0, @isnumeric);
addParameter(p, 'BIN_WIDTH', 50, @isnumeric);
addParameter(p, 'BIN_STEP', 10, @isnumeric);
addParameter(p, 'N_BINS', 24, @isnumeric);
addParameter(p, 'ALPHA', 0.05, @isnumeric);
addParameter(p, 'MIN_RSQ', 0.2, @isnumeric);
addParameter(p, 'MAX_AMP', 500, @isnumeric);
addParameter(p, 'ANGLE', NaN, @isnumeric);
addParameter(p, 'FIG_PATH', [], @ischar);
addParameter(p, 'SAVE_PDF', false, @islogical);

parse(p, data, varargin{:});
data = p.Results.data;
PROBES = p.Results.PROBES;
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
ANGLE = p.Results.ANGLE;
FIG_PATH = p.Results.FIG_PATH;
SAVE_PDF = p.Results.SAVE_PDF;

fprintf('\n------------------------------\n')
if ischar(data)
    [~, filename, ~] = fileparts(data);
    load(data,'S');
    fprintf(sprintf('\n----Data loaded for %s----\n',filename))
else
    S = data;
end

if ~isfield(S, 'sorting')
    error('S.sorting not found -- spike sorting data is required to compute firing rates.');
end
if isempty(PROBES)
    PROBES = unique([S.sorting.probe_index]);
end

[thisTbl, alignStr, ~, ~] = prep_rfsaTable(S, ALIGN_TO, ALIGN_IND, TIME_BIN, ANGLE);

popResults = table();
for pi = 1:numel(PROBES)
    PROBE_INDEX = PROBES(pi);

    prb_name = sprintf('spiketimes_%d', PROBE_INDEX);
    if ~ismember(prb_name, thisTbl.Properties.VariableNames)
        error('Column ''%s'' not found -- these rfsa tables don''t have spike times merged in for PROBE_INDEX %d.', prb_name, PROBE_INDEX);
    end

    sortEntry = S.sorting([S.sorting.probe_index] == PROBE_INDEX);
    clusts = sortEntry.clusters.cluster_id;
    if ismember('best_channel', sortEntry.clusters.Properties.VariableNames)
        chans = sortEntry.clusters.best_channel;
    else
        chans = nan(numel(clusts),1);
    end

    fprintf('\n--- PROBE %d: analyzing %d clusters ---\n', PROBE_INDEX, numel(clusts));
    ticStart = tic;
    for u = 1:numel(clusts)
        clust = clusts(u);

        [frsRaw, bin_edges, xvals, yvals] = format_tableToRFMap(thisTbl, 'PROBE_INDEX', PROBE_INDEX, ...
            'UNITS', (clust+1), 'FIRST_BIN', FIRST_BIN, 'BIN_WIDTH', BIN_WIDTH, 'BIN_STEP', BIN_STEP, 'N_BINS', N_BINS);
        frsRaw = frsRaw{1};
        meanFRs = cellfun(@mean, frsRaw);

        % Excitation and inhibition each get their own independently-detected
        % peak latency bin (largest deviation from that bin's own grid-wide
        % mean, in that direction) -- see ia_rfmpStim_RFMapFit for details
        gridMeanPerBin = squeeze(mean(mean(meanFRs,1),2));
        devFromMean = meanFRs - reshape(gridMeanPerBin,1,1,[]);
        [~, excBinIdx] = max(squeeze(max(max(devFromMean,[],1),[],2)));
        [~, inhBinIdx] = min(squeeze(min(min(devFromMean,[],1),[],2)));

        [hasExcRF, excParams, excRsq, excP] = analyzeRFDirection(meanFRs, xvals, yvals, excBinIdx, 'exc', ALPHA, MIN_RSQ, MAX_AMP);
        [hasInhRF, inhParams, inhRsq, inhP] = analyzeRFDirection(meanFRs, xvals, yvals, inhBinIdx, 'inh', ALPHA, MIN_RSQ, MAX_AMP);

        popResults = [popResults; table(PROBE_INDEX, clust, chans(u), ...
            excBinIdx, bin_edges{excBinIdx}(1), bin_edges{excBinIdx}(2), excP, hasExcRF, excParams(1), excParams(2), excParams(3), excParams(4), excParams(5), excRsq, ...
            inhBinIdx, bin_edges{inhBinIdx}(1), bin_edges{inhBinIdx}(2), inhP, hasInhRF, inhParams(1), inhParams(2), inhParams(3), inhParams(4), inhParams(5), inhRsq, ...
            'VariableNames', {'probe','clust','chan', ...
            'excBinIdx','excLatencyStart','excLatencyEnd','excPValue','hasExcRF','excAmp','excX0','excY0','excSigX','excSigY','excRsq', ...
            'inhBinIdx','inhLatencyStart','inhLatencyEnd','inhPValue','hasInhRF','inhAmp','inhX0','inhY0','inhSigX','inhSigY','inhRsq'})]; %#ok<AGROW>

        if mod(u, 100) == 0
            fprintf('  probe %d: %d / %d clusters done (%.1f min elapsed)\n', PROBE_INDEX, u, numel(clusts), toc(ticStart)/60);
        end
    end
    fprintf('--- PROBE %d complete: %d clusters, %.1f min ---\n', PROBE_INDEX, numel(clusts), toc(ticStart)/60);
end
fprintf('\n------------------------------\n')

%% Population summary figure
hwhmFactor = sqrt(2*log(2));
theta = linspace(0, 2*pi, 100);
latEdges = FIRST_BIN : BIN_STEP : (FIRST_BIN + (N_BINS-1)*BIN_STEP + BIN_WIDTH);

fig = figure('Name','RF Population Summary','Color','w');
fig.Position = [50 50 400*numel(PROBES)+200 900];
tl = tiledlayout(2, numel(PROBES), 'TileSpacing','compact','Padding','compact');

for pi = 1:numel(PROBES)
    PROBE_INDEX = PROBES(pi);
    probeRows = popResults(popResults.probe==PROBE_INDEX, :);
    nClusters = sum(probeRows.hasExcRF | probeRows.hasInhRF);

    % Top row: every detected RF's half-max contour, overlaid
    nexttile(pi);
    hold on;
    for r = 1:height(probeRows)
        if probeRows.hasExcRF(r)
            plot(probeRows.excX0(r) + hwhmFactor*probeRows.excSigX(r)*cos(theta), ...
                 probeRows.excY0(r) + hwhmFactor*probeRows.excSigY(r)*sin(theta), 'r-');
        end
        if probeRows.hasInhRF(r)
            plot(probeRows.inhX0(r) + hwhmFactor*probeRows.inhSigX(r)*cos(theta), ...
                 probeRows.inhY0(r) + hwhmFactor*probeRows.inhSigY(r)*sin(theta), 'b-');
        end
    end
    hold off;
    axis square; grid on;
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    title(sprintf('Probe %d RF Positions (N = %d clusters)', PROBE_INDEX, nClusters));
    prettyFig;

    % Bottom row: latency distributions
    nexttile(numel(PROBES)+pi);
    hold on;
    histogram(probeRows.excLatencyStart(probeRows.hasExcRF), latEdges, 'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none');
    histogram(probeRows.inhLatencyStart(probeRows.hasInhRF), latEdges, 'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none');
    hold off;
    xlabel('Latency (ms)'); ylabel('Count');
    legend({'Excitatory','Inhibitory'}, 'Location','best');
    title(sprintf('Probe %d RF Latencies', PROBE_INDEX));
    prettyFig;
end

title(tl, sprintf('%s --- RF Population Summary (TIME_BIN = [%d, %d] ms rel. to %s)', ...
    S.sess_name, TIME_BIN(1), TIME_BIN(2), alignStr), 'Interpreter', 'none', 'FontSize', 16);

if ~isempty(FIG_PATH)
    if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end
    fileBase = sprintf('RFPopulationSummary_%s_%dto%dms', ALIGN_TO, TIME_BIN(1), TIME_BIN(2));
    if SAVE_PDF
        savebigPDF(fig, fullfile(FIG_PATH, [fileBase '.pdf']));
    else
        savebigPNG(fig, fullfile(FIG_PATH, [fileBase '.png']));
    end
end

end
