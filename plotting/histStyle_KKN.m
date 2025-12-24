function [datamean, datastd] = histStyle_KKN(values,varargin)
    % histStyle_KKN does blah blah
    %
    % Usage:
    %   [datamean, datastd] = histStyle_KKN(values)     
    %   [datamean, datastd] = histStyle_KKN(values, 'FIG_NAME', fig_name)   
    %
    % Inputs:
    %   values - 
    %   'BIN_WIDTH' - optional 
    %   'FIG_TITLE' - optional 
    %   'X_LABEL' - optional 
    %   'Y_LABEL' - optional 
    %   'X_LIMITS' - optional 
    %   'Y_LIMITS' - optional 
    %   'LINE_COLOR' - optional 
    %   'FACE_COLOR' - optional 
    %
    % Output:
    %   rgb_color - Nx3 matrix of RGB values in [0, 1] for each input direction

    p = inputParser;
    addRequired(p, 'values', @isnumeric);
    addParameter(p, 'BIN_WIDTH', [], @isnumeric);
    addParameter(p, 'FIG_NAME', [], @ischar);
    addParameter(p, 'X_LABEL', [], @ischar);
    addParameter(p, 'Y_LABEL', [], @ischar);
    addParameter(p, 'X_LIMITS', [],  @(x) (isnumeric(x)) && numel(x)==2);
    addParameter(p, 'Y_LIMITS', [],  @(x) (isnumeric(x)) && numel(x)==2);
    addParameter(p, 'LINE_COLOR', [],  @(x) (isnumeric(x)) && numel(x)==3);
    addParameter(p, 'FACE_COLOR', [],  @(x) (isnumeric(x)) && numel(x)==3);
    addParameter(p, 'FACE_ALPHA', 1, @isnumeric);
    addParameter(p, 'TEXT_X_POS', [], @isnumeric);
    addParameter(p, 'TEXT_COLOR', [],  @(x) (isnumeric(x)) && numel(x)==3);
    
    parse(p, values, varargin{:});
    values = p.Results.values;
    BIN_WIDTH = p.Results.BIN_WIDTH;
    FIG_NAME = p.Results.FIG_NAME;
    X_LABEL = p.Results.X_LABEL;
    Y_LABEL = p.Results.Y_LABEL;
    X_LIMITS = p.Results.X_LIMITS;
    Y_LIMITS = p.Results.Y_LIMITS;
    LINE_COLOR = p.Results.LINE_COLOR;
    FACE_COLOR = p.Results.FACE_COLOR;
    FACE_ALPHA = p.Results.FACE_ALPHA;
    TEXT_X_POS = p.Results.TEXT_X_POS;
    TEXT_COLOR = p.Results.TEXT_COLOR;

    if isempty(LINE_COLOR)
        LINE_COLOR = [0,0,0];
    end

    if isempty(TEXT_COLOR)
        TEXT_COLOR = [0,0,0];
    end

    if isempty(FACE_COLOR)
        FACE_COLOR = [255,255,255]./255;
    end

    if isempty(BIN_WIDTH)
        val = (abs(max(values)-min(values)))/10;
        BIN_WIDTH = round(val .* (val <= 1) * 10) / 10 + round(val) .* (val > 1);
    end

    % Plot histogram with white bars, range -1 to 
    h = histogram(values, 'binwidth', BIN_WIDTH, 'linewidth', 2, 'facecolor', FACE_COLOR, 'FaceAlpha', FACE_ALPHA);
    hold on;

    % Plot star at MEAN, 0.5 above max of ylim
    plot (mean(values,'omitnan'), max(h.Values)+max(h.Values)/15, '*', 'markersize', 17, 'color', 'k')
    line([mean(values,'omitnan'), mean(values,'omitnan') ], [0 max(h.Values)+max(h.Values)./15], 'color', LINE_COLOR, 'linewidth', 3)

    if ~isempty(X_LIMITS)
        xlim(X_LIMITS)
    end
    if ~isempty(Y_LIMITS)
        ylim(Y_LIMITS)
    else
        Y_LIMITS = [0 max(h.Values)];
    end

    axis square
    
    % draw thin vertical dotted line at 0 for reference
    if h.BinLimits(1)<0 && h.BinLimits(2)>0
        line([0,0], Y_LIMITS + Y_LIMITS./10, 'color', 'k', 'linewidth', 2, 'linestyle', ':')
    end
    
    % draw thick vertical solid line at mean (redundant with star above)
    
    if ~isempty(FIG_NAME)
        title(FIG_NAME,'Interpreter','tex') 
    end

    if ~isempty(X_LABEL)
        xlabel(X_LABEL) 
    end
    if ~isempty(Y_LABEL)
        ylabel(Y_LABEL) 
    end              
    
    prettyFig;
    
    % Perform signrank test to determine if distribution is significantly
    % different from 0 (Note the technical wrinkle that the signrank test used
    % the median but we are displaying the mean; this is standard)

    datamean = mean(values,'omitnan');
    datastd = std(values,'omitnan');
    [pvalue, hypothesis] = signrank(values);

    if isempty(TEXT_X_POS)
        TEXT_X_POS = h.BinLimits(2) - (h.BinLimits(2)-h.BinLimits(1))/5;
    end
    
    % Plot values on figure
    if hypothesis == 1 % if null hypothesis is rejected/ statistically signif, use green font for pvalue
        text (TEXT_X_POS, max(Y_LIMITS)*0.97, ['n = ', num2str(length(values))], 'fontsize', 15, 'color', TEXT_COLOR ) % sample size
        text (TEXT_X_POS, max(Y_LIMITS)*0.91, ['mean = ', sprintf('%0.3f', datamean)],'fontsize', 15, 'color', TEXT_COLOR ) % mean
        text (TEXT_X_POS, max(Y_LIMITS)*0.85, ['std = ', sprintf('%0.3f', datastd)], 'fontsize', 15, 'color', TEXT_COLOR ) % std
        text (TEXT_X_POS, max(Y_LIMITS)*0.79, ['p = ', sprintf('%0.3f', pvalue)], 'fontweight', 'bold', 'fontsize', 15, 'color', TEXT_COLOR, 'fontweight', 'bold' ) % pvalue
    else % If not significant, pvalue in black
        text (TEXT_X_POS, max(Y_LIMITS)*0.97, ['n = ', num2str(length(values))], 'fontsize', 15, 'color', TEXT_COLOR ) % sample size
        text (TEXT_X_POS, max(Y_LIMITS)*0.91, ['mean = ', sprintf('%0.3f', datamean)], 'fontsize', 15, 'color', TEXT_COLOR ) % mean
        text (TEXT_X_POS, max(Y_LIMITS)*0.85, ['std = ', sprintf('%0.3f', datastd)], 'fontsize', 15, 'color', TEXT_COLOR ) % std
        text (TEXT_X_POS, max(Y_LIMITS)*0.79, ['p = ', sprintf('%0.3f', pvalue)], 'fontsize', 15, 'color', TEXT_COLOR ) % pvalue
    end
end

