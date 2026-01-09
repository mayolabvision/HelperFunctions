function [bin_idx, bin_n, bin_centers, bin_edges] = binXvalues(x, varargin)
% BINXVALUES Bin a numeric vector into equal-width or quantile bins
%
% [bin_idx, bin_n, bin_centers, bin_edges] = binXvalues(x, 'Name', Value, ...)
%
% Inputs:
%   x         - numeric vector (Nx1)
% Name-Value pairs:
%   'numBins' - number of bins (default: 8)
%   'BinType' - 'equalwidth' (default) or 'quantile'
%
% Outputs:
%   bin_idx     - Nx1 array, bin index for each value in x
%   bin_n       - numBins x 1 array, number of values in each bin
%   bin_centers - numBins x 1 array, center value of each bin
%   bin_edges   - edges of the bins (numBins+1 x 1)

% Parse inputs
p = inputParser;
addParameter(p,'NumBins',8,@(v) isnumeric(v) && isscalar(v) && v>0);
addParameter(p,'BinType','equalwidth',@(s) any(validatestring(s,{'equalwidth','quantile'})));
parse(p,varargin{:});
numBins = p.Results.NumBins;
binType = p.Results.BinType;

x = x(:);  % ensure column

switch binType
    case 'equalwidth'
        bin_edges = linspace(min(x), max(x), numBins+1);
    case 'quantile'
        bin_edges = quantile(x, linspace(0,1,numBins+1));
end

% Assign each value to a bin
bin_idx = discretize(x, bin_edges);

% Compute bin statistics
bin_n = accumarray(bin_idx, 1, [numBins,1]);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

end