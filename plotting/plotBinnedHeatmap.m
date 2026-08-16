function [n_bins, bin_extrema, heatMat] = plotBinnedHeatmap(x, y, c, varargin)
% PLOTBINNEDHEATMAP Plots a heatmap from raw x, y values binned by edges
%
% n_bins = plotBinnedHeatmap(x, y, c, 'Name', Value, ...)

% Parse inputs
p = inputParser;
addParameter(p,'xEdges',[]);
addParameter(p,'yEdges',[]);
addParameter(p,'nBins', 8, @isnumeric);
addParameter(p,'Type','mean',@(s) any(validatestring(s,{'mean','median','std','count'})));
addParameter(p,'cmap',parula);
addParameter(p,'climit',[]);
addParameter(p,'scaleBins',true,@islogical);
addParameter(p, 'RETURN_CELL',false,@islogical);

parse(p,varargin{:});
xEdges = p.Results.xEdges;
yEdges = p.Results.yEdges;
nBins = p.Results.nBins;
aggType = p.Results.Type;
cmap = p.Results.cmap;
climit = p.Results.climit;
scaleBins = p.Results.scaleBins;
RETURN_CELL = p.Results.RETURN_CELL;

x = x(:); y = y(:);
if nargin < 3 || isempty(c)
    c = ones(size(x)); % default to count
    aggType = 'count';
end
c = c(:);

% Define default bin edges if not provided
if isempty(xEdges)
    xEdges = linspace(min(x), max(x), nBins); % 8 bins
end
if isempty(yEdges)
    yEdges = linspace(min(y), max(y), nBins); 
end

% Bin indices for each point
[~,~,xBin] = histcounts(x, xEdges);
[~,~,yBin] = histcounts(y, yEdges);

% Compute aggregated value per bin
if RETURN_CELL
    heatMat = cell(length(yEdges)-1, length(xEdges)-1);
    n_bins  = zeros(length(yEdges)-1, length(xEdges)-1);
    for ix = 1:length(xEdges)-1
        for iy = 1:length(yEdges)-1
            inBin = xBin == ix & yBin == iy;
            n_bins(iy, ix) = sum(inBin); % store count
            if sum(inBin) > 0
                heatMat{iy, ix} = c(inBin); % store values
            end
        end
    end

    bin_extrema = [];
    return
end

heatMat = nan(length(yEdges)-1, length(xEdges)-1); % initialize as NaN
n_bins  = zeros(length(yEdges)-1, length(xEdges)-1);

for ix = 1:length(xEdges)-1
    for iy = 1:length(yEdges)-1
        inBin = xBin == ix & yBin == iy;
        vals = c(inBin);
        n_bins(iy, ix) = sum(inBin); % store count

        if isempty(vals) % <-- no points in bin
            heatMat(iy, ix) = NaN; % will appear white
        else
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
end

% ---- Find min / max bins (ignore NaNs) ----
validVals = heatMat(~isnan(heatMat));

[minVal, ~] = min(validVals);
[maxVal, ~] = max(validVals);

[minIy, minIx] = find(heatMat == minVal, 1);
[maxIy, maxIx] = find(heatMat == maxVal, 1);

bin_extrema.min.value = minVal;
bin_extrema.min.bin_index = [minIy, minIx];
bin_extrema.min.x_edges = xEdges(minIx:minIx+1);
bin_extrema.min.y_edges = yEdges(minIy:minIy+1);

bin_extrema.max.value = maxVal;
bin_extrema.max.bin_index = [maxIy, maxIx];
bin_extrema.max.x_edges = xEdges(maxIx:maxIx+1);
bin_extrema.max.y_edges = yEdges(maxIy:maxIy+1);

% Plot heatmap
% ----- Plot heatmap -----
if scaleBins
    % TRUE geometric bin scaling
    [Xgrid, Ygrid] = meshgrid(xEdges, yEdges);

    % pcolor requires one extra row/column
    heatMat_plot = padarray(heatMat,[1 1],NaN,'post');

    h = pcolor(Xgrid, Ygrid, heatMat_plot);
    shading flat
else
    % Original equal-sized bins
    h = imagesc(xEdges, yEdges, heatMat);
    set(gca,'YDir','normal');
end

% handle NaNs as white
colormap(cmap); % prepend white to colormap
clim([minVal maxVal])   % keep color limits consistent
colorbar;

% Apply color limits if provided
if ~isempty(climit)
    clim(climit);
end

end