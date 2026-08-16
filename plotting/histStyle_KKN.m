function [datamean, datastd] = histStyle_KKN(values, varargin)
% histStyle_KKN  Plot a styled histogram with summary statistics overlay.
%
% Usage:
%   [datamean, datastd] = histStyle_KKN(values)
%   [datamean, datastd] = histStyle_KKN(values, 'Name', value, ...)
%
% Required input:
%   values       - numeric vector
%
% Optional name-value parameters:
%   'BIN_WIDTH'  - scalar bin width (default: auto, ~10 bins); ignored when LOG_X is true
%   'NORM'       - normalization passed to histogram() (default: 'count')
%   'FIG_NAME'   - axes title string
%   'X_LABEL'    - x-axis label
%   'Y_LABEL'    - y-axis label
%   'X_LIMITS'   - 1×2 x-axis limits
%   'Y_LIMITS'   - 1×2 y-axis limits
%   'LINE_COLOR' - 1×3 RGB for mean line and marker (default: black)
%   'FACE_COLOR' - 1×3 RGB bar fill (default: white)
%   'FACE_ALPHA' - bar transparency 0–1 (default: 1)
%   'TEXT_X_POS' - x position of stats text block
%   'TEXT_Y_POS' - y position of stats text block
%   'TEXT_COLOR' - 1×3 RGB for stats text (default: black)
%   'LOG_X'      - logical, use log-spaced bins and log x-axis (default: false)  % NEW
%
% Outputs:
%   datamean - mean of values, NaN omitted
%   datastd  - std of values, NaN omitted

    p = inputParser;
    addRequired(p,  'values',     @isnumeric);
    addParameter(p, 'BIN_WIDTH',  [],        @isnumeric);
    addParameter(p, 'NORM',       'count',   @ischar);
    addParameter(p, 'FIG_NAME',   [],        @ischar);
    addParameter(p, 'X_LABEL',    [],        @ischar);
    addParameter(p, 'Y_LABEL',    [],        @ischar);
    addParameter(p, 'X_LIMITS',   [],        @(x) isnumeric(x) && numel(x)==2);
    addParameter(p, 'Y_LIMITS',   [],        @(x) isnumeric(x) && numel(x)==2);
    addParameter(p, 'LINE_COLOR', [],        @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'FACE_COLOR', [],        @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'FACE_ALPHA', 1,         @isnumeric);
    addParameter(p, 'TEXT_X_POS', [],        @isnumeric);
    addParameter(p, 'TEXT_Y_POS', [],        @isnumeric);
    addParameter(p, 'TEXT_COLOR', [],        @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'LOG_X',      false,     @islogical);    % NEW
    parse(p, values, varargin{:});

    values     = p.Results.values;
    BIN_WIDTH  = p.Results.BIN_WIDTH;
    NORM       = p.Results.NORM;
    FIG_NAME   = p.Results.FIG_NAME;
    X_LABEL    = p.Results.X_LABEL;
    Y_LABEL    = p.Results.Y_LABEL;
    X_LIMITS   = p.Results.X_LIMITS;
    Y_LIMITS   = p.Results.Y_LIMITS;
    LINE_COLOR = p.Results.LINE_COLOR;
    FACE_COLOR = p.Results.FACE_COLOR;
    FACE_ALPHA = p.Results.FACE_ALPHA;
    TEXT_X_POS = p.Results.TEXT_X_POS;
    TEXT_Y_POS = p.Results.TEXT_Y_POS;
    TEXT_COLOR = p.Results.TEXT_COLOR;
    LOG_X      = p.Results.LOG_X;              % NEW

    if isempty(LINE_COLOR), LINE_COLOR = [0 0 0]; end
    if isempty(TEXT_COLOR), TEXT_COLOR = [0 0 0]; end
    if isempty(FACE_COLOR), FACE_COLOR = [1 1 1]; end

    % NEW: log mode uses 20 log-spaced bin edges; BIN_WIDTH is ignored.
    % Values <= 0 are excluded from log plots (they cannot be represented).
    if LOG_X
        pos_vals  = values(values > 0);
        log_edges = logspace(log10(min(pos_vals)), log10(max(pos_vals)), 21);
        h = histogram(pos_vals, log_edges, 'LineWidth', 2, ...
            'FaceColor', FACE_COLOR, 'FaceAlpha', FACE_ALPHA, 'Normalization', NORM);
        set(gca, 'XScale', 'log');
    else
        if isempty(BIN_WIDTH)
            val = (max(values) - min(values)) / 10;
            if val <= 1
                BIN_WIDTH = round(val, 1);   % nearest 0.1
            else
                BIN_WIDTH = round(val);      % nearest integer
            end
            BIN_WIDTH = max(BIN_WIDTH, eps); % guard against zero-range data
        end
        h = histogram(values, 'BinWidth', BIN_WIDTH, 'LineWidth', 2, ...
            'FaceColor', FACE_COLOR, 'FaceAlpha', FACE_ALPHA, 'Normalization', NORM);
    end
    hold on;

    if ~isempty(X_LIMITS), xlim(X_LIMITS); end

    if ~isempty(Y_LIMITS)
        ylim(Y_LIMITS);
    else
        Y_LIMITS = [0, max(h.Values) + max(h.Values)/10];
        ylim(Y_LIMITS);
    end

    datamean = mean(values, 'omitnan');
    datastd  = std(values,  'omitnan');

    % Mean marker and vertical line
    line([datamean, datamean], [0, max(Y_LIMITS) - max(Y_LIMITS)/20], ...
        'Color', LINE_COLOR, 'LineWidth', 3);
    plot(datamean, max(Y_LIMITS) - max(Y_LIMITS)/20, ...
        '*', 'MarkerSize', 17, 'Color', LINE_COLOR);

    axis square;

    % Zero-reference line (only when data spans zero; skipped in log mode
    % since log axes cannot represent zero)
    if ~LOG_X && h.BinLimits(1) < 0 && h.BinLimits(2) > 0
        line([0 0], Y_LIMITS + Y_LIMITS./10, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');
    end

    if ~isempty(FIG_NAME), title(FIG_NAME, 'Interpreter', 'tex'); end
    if ~isempty(X_LABEL),  xlabel(X_LABEL); end
    if ~isempty(Y_LABEL),  ylabel(Y_LABEL); end

    % Stats text overlay. signrank tests whether the distribution is centered
    % at zero — bold p-value when significant.
    [pvalue, hypothesis] = signrank(values);

    % NEW: in log mode, compute text x-position in log space so it lands at
    % the same relative position as it would on a linear axis.
    if isempty(TEXT_X_POS)
        if LOG_X
            TEXT_X_POS = 10^(log10(h.BinLimits(2)) - ...
                (log10(h.BinLimits(2)) - log10(h.BinLimits(1))) / 5);
        else
            TEXT_X_POS = h.BinLimits(2) - (h.BinLimits(2) - h.BinLimits(1)) / 5;
        end
    end
    if isempty(TEXT_Y_POS)
        TEXT_Y_POS = max(Y_LIMITS);
    end

    pval_weight = 'normal';
    if hypothesis == 1, pval_weight = 'bold'; end

    text(TEXT_X_POS, TEXT_Y_POS*0.97, ['n = ',    num2str(length(values))],        'FontSize', 15, 'Color', TEXT_COLOR);
    text(TEXT_X_POS, TEXT_Y_POS*0.91, ['mean = ', sprintf('%0.3f', datamean)],      'FontSize', 15, 'Color', TEXT_COLOR);
    text(TEXT_X_POS, TEXT_Y_POS*0.85, ['std = ',  sprintf('%0.3f', datastd)],       'FontSize', 15, 'Color', TEXT_COLOR);
    text(TEXT_X_POS, TEXT_Y_POS*0.79, ['p = ',    sprintf('%0.3f', pvalue)], ...
        'FontSize', 15, 'Color', TEXT_COLOR, 'FontWeight', pval_weight);
end