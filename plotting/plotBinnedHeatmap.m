function n_bins = plotBinnedHeatmap(x, y, c, varargin)
% PLOTBINNEDHEATMAP Plots a heatmap from raw x, y values binned by edges
%
% n_bins = plotBinnedHeatmap(x, y, c, 'Name', Value, ...)
%
% Inputs:
%   x - Nx1 numeric vector (raw x values)
%   y - Nx1 numeric vector (raw y values)
%   c - Nx1 numeric vector (values to aggregate; optional if only counting)
%
% Name-Value pairs:
%   'xEdges'   - vector of edges for x-axis bins (default: 8 evenly spaced bins)
%   'yEdges'   - vector of edges for y-axis bins (default: 8 evenly spaced bins)
%   'Type'     - 'mean' (default), 'median', 'std', 'count'
%   'cmap'     - colormap for the heatmap (default: parula)
%   'climit'   - 2-element vector to limit the colorbar (default: automatic)

% Parse inputs
p = inputParser;
addParameter(p,'xEdges',[]);
addParameter(p,'yEdges',[]);
addParameter(p,'Type','mean',@(s) any(validatestring(s,{'mean','median','std','count'})));
addParameter(p,'cmap',parula);
addParameter(p,'climit',[]);
parse(p,varargin{:});
xEdges = p.Results.xEdges;
yEdges = p.Results.yEdges;
aggType = p.Results.Type;
cmap = p.Results.cmap;
climit = p.Results.climit;

x = x(:); y = y(:);
if nargin < 3 || isempty(c)
    c = ones(size(x)); % default to count
    aggType = 'count';
end
c = c(:);

% Define default bin edges if not provided
if isempty(xEdges)
    xEdges = linspace(min(x), max(x), 9); % 8 bins
end
if isempty(yEdges)
    yEdges = linspace(min(y), max(y), 9); 
end

% Bin indices for each point
[~,~,xBin] = histcounts(x, xEdges);
[~,~,yBin] = histcounts(y, yEdges);

% Preallocate matrices
heatMat = nan(length(yEdges)-1, length(xEdges)-1);
n_bins  = zeros(length(yEdges)-1, length(xEdges)-1);

% Compute aggregated value per bin
for ix = 1:length(xEdges)-1
    for iy = 1:length(yEdges)-1
        inBin = xBin == ix & yBin == iy;
        vals = c(inBin);
        n_bins(iy, ix) = sum(inBin); % store count

        switch aggType
            case 'mean'
                heatMat(iy, ix) = mean(vals,'omitnan');
            case 'median'
                heatMat(iy, ix) = median(vals,'omitnan');
            case 'std'
                heatMat(iy, ix) = std(vals,'omitnan');
            case 'count'
                heatMat(iy, ix) = sum(inBin);
        end
    end
end

% Plot heatmap
h = imagesc(xEdges, yEdges, heatMat);
set(gca,'YDir','normal'); % make y increase upwards
colormap(cmap);
colorbar;

% Apply color limits if provided
if ~isempty(climit)
    clim(climit);
end

end