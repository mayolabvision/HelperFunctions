function [datamean, datastd] = polarHistStyle_KKN(values,varargin)
    % polarHistStyle_KKN plots a polar histogram of angular data
    %
    % Usage:
    %   [datamean, datastd] = polarHistStyle_KKN(values)     
    %   [datamean, datastd] = polarHistStyle_KKN(values, 'FIG_NAME', fig_name)   
    %
    % Inputs:
    %   values - angular values in DEGREES
    %   'BIN_WIDTH' - optional (bin width in DEGREES)
    %   'FIG_NAME' - optional title
    %   'LINE_COLOR' - optional RGB [0 1] color
    %   'FACE_COLOR' - optional RGB [0 1] color
    %
    % Outputs:
    %   datamean - circular mean (DEGREES)
    %   datastd  - circular std deviation (DEGREES)

    p = inputParser;
    addRequired(p, 'values', @isnumeric);
    addParameter(p, 'BIN_WIDTH', [], @isnumeric);
    addParameter(p, 'FIG_NAME', [], @ischar);
    addParameter(p, 'LINE_COLOR', [],  @(x) (isnumeric(x)) && numel(x)==3);
    addParameter(p, 'FACE_COLOR', [],  @(x) (isnumeric(x)) && numel(x)==3);
    addParameter(p, 'FACE_ALPHA', 1, @isnumeric);
    addParameter(p, 'MARK_RAD', [], @isnumeric);
    addParameter(p, 'RLIM', [], @isnumeric);
    addParameter(p, 'TEXT_Y_POS', [], @isnumeric);
    addParameter(p, 'TEXT_COLOR', [],  @(x) (isnumeric(x)) && numel(x)==3);
    parse(p, values, varargin{:});

    values_deg = p.Results.values(:); % ensure column
    BIN_WIDTH_deg = p.Results.BIN_WIDTH;
    FIG_NAME  = p.Results.FIG_NAME;
    LINE_COLOR = p.Results.LINE_COLOR;
    FACE_COLOR = p.Results.FACE_COLOR;
    MARK_RAD = p.Results.MARK_RAD;
    FACE_ALPHA = p.Results.FACE_ALPHA;
    RLIM = p.Results.RLIM;
    TEXT_Y_POS = p.Results.TEXT_Y_POS;
    TEXT_COLOR = p.Results.TEXT_COLOR;

    if isempty(LINE_COLOR), LINE_COLOR = [0 0 0]; end
    if isempty(FACE_COLOR), FACE_COLOR = [0.8 0.8 0.8]; end
    if isempty(BIN_WIDTH_deg),  BIN_WIDTH_deg = 10; end % default 10°

    % Convert to radians for plotting/stats
    values = deg2rad(values_deg);
    BIN_WIDTH = deg2rad(BIN_WIDTH_deg);

    % --- Plot polar histogram ---
    h = polarhistogram(values, 'BinWidth', BIN_WIDTH, ...
        'FaceColor', FACE_COLOR, 'EdgeColor', LINE_COLOR, 'LineWidth', 1.5, 'FaceAlpha', FACE_ALPHA);
    hold on;

    if ~isempty(RLIM)
        rlim([0 RLIM])
    end

    % --- Circular mean and std (in radians first) ---
    datamean_rad = circ_mean(values);
    datastd_rad  = circ_std(values);

    [pval, ~] = circ_rtest(values);

    % Convert back to degrees for output
    datamean = rad2deg(datamean_rad);
    datastd  = rad2deg(datastd_rad);

    % plot mean direction as a radial line
    if isempty(MARK_RAD)
        rmax = max(h.Values) + max(h.Values)/10;
    else
        rmax = MARK_RAD;
    end

    % Scale line length: variance=0 → rmax, variance=std_max → 0
    line_length = rmax * max(0, 1 - datastd / 360);
    
    polarplot([datamean_rad datamean_rad], [0 line_length], ...
        'Color', LINE_COLOR, 'LineWidth', 2);

    % mark mean direction with a star
    polarplot(datamean_rad,line_length,'p','MarkerSize',12,'MarkerFaceColor',LINE_COLOR,'MarkerEdgeColor',LINE_COLOR)

    % optional title
    if ~isempty(FIG_NAME)
        title(FIG_NAME)
    end

    if isempty(TEXT_Y_POS)
        TEXT_Y_POS = 1.05;
    end

    if isempty(TEXT_COLOR)
        TEXT_COLOR = [0,0,0];
    end

    % --- Add stats text in top-right corner (normalized coords) ---
    ax = gca;
    text(ax, 1, TEXT_Y_POS+0.05, sprintf('n = %d', numel(values_deg)), ...
        'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color',TEXT_COLOR);
    text(ax, 1, TEXT_Y_POS, sprintf('mean = %.2f°', wrapTo360(datamean)), ...
        'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color',TEXT_COLOR);
    text(ax, 1, TEXT_Y_POS-0.05, sprintf('std = %.2f°', datastd), ...
        'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color',TEXT_COLOR);

    if pval < 0.05
        text(ax, 1, TEXT_Y_POS-0.1, sprintf('pval = %.3f°', pval), ...
        'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color',TEXT_COLOR,'FontWeight','bold');
    else
        text(ax, 1, TEXT_Y_POS-0.1, sprintf('p = %.3f°', pval), ...
        'Units','normalized','HorizontalAlignment','left','FontSize',12,'Color',TEXT_COLOR,'FontWeight','normal');
    end

    % summary stats printed in Command Window too
end