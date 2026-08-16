function [bin_mean, bin_sem, bin_n, bin_centers, bin_edges] = binYbyX(x, y, varargin)
% BINYBYX Bin x-values and compute mean/SEM of y in each bin
%
% Inputs
%   x, y        : column vectors (same length)
%
% Optional name-value pairs
%   'NumBins'   : number of bins (default = 10)
%   'BinType'   : 'equalwidth' or 'quantile' (default = 'equalwidth')
%   'BinEdges'  : explicit bin edges (vector). Overrides NumBins & BinType.
%
% Outputs
%   bin_mean    : mean y per bin
%   bin_sem     : SEM of y per bin
%   bin_n       : number of points per bin
%   bin_centers : center x-value of each bin
%   bin_edges   : bin edges used

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('x', @(v)isvector(v) && isnumeric(v));
p.addRequired('y', @(v)isvector(v) && isnumeric(v));
p.addParameter('NumBins', 10, @(v)isnumeric(v) && isscalar(v) && v>0);
p.addParameter('BinType', 'equalwidth', @(s)ischar(s) || isstring(s));
p.addParameter('BinEdges', [], @(v)isnumeric(v) && isvector(v));
p.parse(x, y, varargin{:});

nBins    = p.Results.NumBins;
binType  = lower(p.Results.BinType);
binEdges = p.Results.BinEdges;

x = x(:);
y = y(:);

assert(numel(x) == numel(y), 'x and y must be the same length');

% remove NaNs
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

% -------------------- define bins --------------------
if ~isempty(binEdges)
    % Explicit edges override everything
    bin_edges = binEdges(:);
    nBins = numel(bin_edges) - 1;

else
    switch binType
        case 'equalwidth'
            bin_edges = linspace(min(x), max(x), nBins+1);

        case 'quantile'
            qs = linspace(0,1,nBins+1);
            bin_edges = quantile(x, qs);
            bin_edges(1)   = min(x);
            bin_edges(end) = max(x);

        otherwise
            error('BinType must be ''equalwidth'' or ''quantile''');
    end
end

% -------------------- bin data --------------------
bin_mean    = nan(nBins,1);
bin_sem     = nan(nBins,1);
bin_centers = nan(nBins,1);
bin_n       = zeros(nBins,1);

for i = 1:nBins
    if i < nBins
        inBin = x >= bin_edges(i) & x < bin_edges(i+1);
    else
        inBin = x >= bin_edges(i) & x <= bin_edges(i+1);
    end

    yi = y(inBin);
    bin_n(i) = numel(yi);

    if bin_n(i) > 0
        bin_mean(i) = mean(yi);
        bin_sem(i)  = std(yi) / sqrt(bin_n(i));
        bin_centers(i) = mean(bin_edges(i:i+1));
    end
end
end