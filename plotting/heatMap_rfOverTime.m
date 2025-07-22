function tl = heatMap_rfOverTime(frs, varargin)
    % Generates a heatmap visualization of receptive fields (RF) over time.
    %
    % This function takes an array of firing rates (FRs) measured over a grid 
    % of spatial bins and time bins and plots a heatmap for each time window.
    % It is useful for analyzing how the RF of a neuron evolves before and 
    % during the presentation of stimuli in space. 
    %
    %%% Required Input: %%%
    %   frs        -  3D data structure representing mean firing rates.
    %                 Can be either:
    %                 - A 3D numeric array (X bins, Y bins, Time bins).
    %                   Example: 12x12x24 (12x12 spatial grid, 24 time bins).
    %                 - A 3D cell array where each cell contains an array 
    %                   of firing rates, which will be averaged to obtain 
    %                   a numeric array of mean firing rates.
    %
    %%% Optional Parameters: %%%
    %   'BIN_EDGES'  -  Cell array (1xN) specifying time bin edges.
    %                   Each cell contains a [start, end] time (e.g., [0 50]).
    %                   Default: not specified.
    %
    %   'X_VALS'     -  Numeric array for x-axis labels (spatial positions).
    %                   Example: [-30 -25 -20 ... 25 30] (degrees).
    %                   Default: 1:size(frs,1).
    %
    %   'Y_VALS'     -  Numeric array for y-axis labels (spatial positions).
    %                   Example: [-30 -25 -20 ... 25 30] (degrees).
    %                   Default: 1:size(frs,2).
    %
    %   'INTERP'     -  Logical flag (true/false) for enabling heatmap interpolation.
    %                   If true, smooths the heatmap using interpolation.
    %                   Default: false.
    %
    %%% Output: %%%
    %   tl  -  Tiled layout object, allowing additional customization
    %          (e.g., adding titles or annotations).
    %
    %%% Example Usage: %%%
    %   tl = heatMap_rfOverTime(frs,'BIN_EDGES',bin_edges, 'INTERP', false,...
    %                           'X_VALS',xvals, 'Y_VALS', yvals);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    defaultXVals = 1:size(frs,1);
    defaultYVals = 1:size(frs,2);
    defaultInterp = false;

    % Create an input parser
    p = inputParser;
    addRequired(p, 'frs', @(x) ((isnumeric(x) || iscell(x)) && numel(size(x)) == 3));
    addParameter(p, 'BIN_EDGES', [], @(x) (iscell(x) && numel(x) == size(frs,3)));
    addParameter(p, 'X_VALS', defaultXVals, @(x) (isnumeric(x) && numel(x) == size(frs,1)));
    addParameter(p, 'Y_VALS', defaultYVals, @(x) (isnumeric(x) && numel(x) == size(frs,2)));
    addParameter(p, 'INTERP', defaultInterp, @(x) (islogical(x)));

    % Parse the inputs
    parse(p, frs, varargin{:});

    % Assign parsed values to variables
    frs = p.Results.frs;
    bin_edges = p.Results.BIN_EDGES;
    xvals= p.Results.X_VALS;
    yvals = p.Results.Y_VALS;
    interp = p.Results.INTERP;

    if iscell(frs)
        frs = cellfun(@(q) mean(q), frs);
    end

    if size(xvals,1)>size(xvals,2)
        xvals = xvals';
        yvals = yvals';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tl = tiledlayout(4,ceil(numel(bin_edges)/4));
    tl.TileSpacing = 'compact';
    tl.Padding = 'tight';

    ax = gobjects(1, length(bin_edges)); 
    for bin = 1:length(bin_edges)
        ax(bin) = nexttile(tl);
        
        % heatmap(xvals,yvals,frs(:,:,bin), ...
        %     'Colormap', parula, ...
        %     'ColorbarVisible', 'on', ...
        %     'XLabel', 'X Position', ...
        %     'YLabel', 'Y Position', ...
        %     'Title', 'Mean Firing Rate');
        

        frs2 = [frs(:,:,bin); frs(end,:,bin)];
        frs2 = [frs2 frs2(:,end)];

        dx = mean(diff(xvals));
        dy = mean(diff(yvals));
        xedges = [xvals - dx/2, xvals(end) + dx/2];
        yedges = [yvals - dy/2, yvals(end) + dy/2];

        %hmap = pcolor(xvals, yvals, frs(:,:,bin));
        hmap = pcolor(xedges, yedges, frs2);

        if interp
            hmap.FaceColor = 'interp';
        end
        
        if ~isempty(bin_edges)
            title(sprintf('[%d , %d) ms',bin_edges{bin}(1),bin_edges{bin}(2)))
        end
    
        axis square;
        prettyFig;
    end
    
    if max(frs(:))>0
        set(ax, 'Colormap', gray, 'CLim', [min(frs(:)) max(frs(:))])
    end
    cbh = colorbar(ax(end)); 
    cbh.Layout.Tile = 'east'; 
    cbh.Label.String = 'Firing Rate (Hz)';
    cbh.Label.Rotation = 270;

    xlabel(tl,'Horizontal (deg)','fontsize',16)
    ylabel(tl,'Vertical (deg)','fontsize',16)

end