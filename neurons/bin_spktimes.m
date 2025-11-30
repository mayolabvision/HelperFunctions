function spks_binned = bin_spktimes(spktimes, varargin)
% bin_spktimes - Bin spike times into fixed-size or sliding bins.
%
% Usage:
%   spks_binned = bin_spktimes(spktimes)
%   spks_binned = bin_spktimes(spktimes, 'StartTime', 0, 'EndTime', 1000, 'BinSize', 50)
%   spks_binned = bin_spktimes(spktimes, 'BinSize', 50, 'SlidingBin', 25)
%
% Defaults:
%   StartTime  = 0
%   EndTime    = max(spktimes)
%   BinSize    = 50 ms   (window width)
%   SlidingBin = BinSize (step size → no overlap)

    % ---- parser ----
    p = inputParser;
    addRequired(p, 'spktimes', @(x) isnumeric(x));

    addParameter(p, 'StartTime', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'EndTime', max(spktimes), @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'BinSize', 50, @(x) isnumeric(x) && isscalar(x));

    % Sliding step (default = BinSize → no overlap)
    addParameter(p, 'SlidingBin', [], @(x) isnumeric(x) && isscalar(x));

    parse(p, spktimes, varargin{:});

    startTime  = p.Results.StartTime;
    endTime    = p.Results.EndTime;
    binSize    = p.Results.BinSize;
    slideSize  = p.Results.SlidingBin;

    % set default sliding step equal to bin width
    if isempty(slideSize)
        slideSize = binSize;
    end

    % ---- sliding window edges ----
    % window starts
    binStarts = startTime : slideSize : (endTime - binSize);

    % preallocate counts
    counts = zeros(1, numel(binStarts));

    % ---- count spikes per window ----
    for i = 1:numel(binStarts)
        bStart = binStarts(i);
        bEnd   = bStart + binSize;
        counts(i) = sum(spktimes >= bStart & spktimes < bEnd);
    end

    spks_binned = counts;
end