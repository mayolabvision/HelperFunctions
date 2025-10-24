function spks_binned = bin_spktimes(spktimes, varargin)
% bin_spktimes - Bin spike times into fixed-size bins.
%
% Usage:
%   spks_binned = bin_spktimes(spktimes)
%   spks_binned = bin_spktimes(spktimes, 'StartTime', 0, 'EndTime', 1000, 'BinSize', 50)
%
% Defaults:
%   StartTime = 0
%   EndTime   = max(spktimes)
%   BinSize   = 50

    % set up parser
    p = inputParser;
    addRequired(p, 'spktimes', @(x) isnumeric(x));
    addParameter(p, 'StartTime', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'EndTime', max(spktimes), @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'BinSize', 50, @(x) isnumeric(x) && isscalar(x));

    parse(p, spktimes, varargin{:});
    startTime = p.Results.StartTime;
    endTime   = p.Results.EndTime;
    binSize   = p.Results.BinSize;

    % binning
    edges = startTime:binSize:endTime;
    [counts, ~] = histcounts(spktimes, edges);
    spks_binned = counts;
end