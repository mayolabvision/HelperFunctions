function sdf = make_sdf(spike_times, tWin, sigma)
% make_sdf — Compute spike density functions by Gaussian smoothing.
%
% INPUTS:
%   spike_times : A cell array where each cell contains a *vector* of spike
%                 times (in ms) for one trial.
%                 • Required type:  N×1 or 1×N cell array
%                 • Each element:   1×M or M×1 numeric vector of spike times
%                 • Spike times may be negative/positive (relative to event)
%                 Example:
%                     spike_times{1} = [ -20, 5, 17, 33, 102 ];
%                     spike_times{2} = [ 7, 9, 11, 200 ];
%
%   tWin        : Two-element vector [tStart tEnd] defining the output time axis
%                 in ms. Example: [-300 500].
%
%   sigma       : Standard deviation (in ms) of the Gaussian smoothing kernel.
%                 Common values = 10 ms.
%
% OUTPUT:
%   sdf         : A cell array of the same size as spike_times, where each cell
%                 contains a 1×T vector (T = tEnd - tStart + 1) representing the
%                 spike density function (in spikes/sec) for that trial.
%
% METHOD:
%   Converts each trial's spike times into a 1-ms resolution binary spike train,
%   then convolves that train with a normalized Gaussian kernel to obtain an SDF.
%
% NOTES:
%   • Sampling resolution is fixed at 1 ms.
%   • Output units are spikes/second.
%   • Convolution is applied using 'same' to preserve the original window size.

    % Sampling rate (1 ms resolution)
    Fs = 1000;

    % Create time vector for the output SDF
    tvec = tWin(1) : tWin(2);

    % ---- Gaussian kernel centered around zero ----
    % Kernel extends ±4 sigma
    kt = -4*sigma : 4*sigma;

    % Unnormalized Gaussian
    k = exp(-(kt.^2) / (2 * sigma^2));

    % Normalize to area = 1 (so convolution yields spikes/ms)
    k = k / sum(k);

    % Preallocate output SDF as a cell array of matching size
    sdf = cell(size(spike_times));

    % ---- Compute SDF for each trial ----
    for i = 1:numel(spike_times)

        % Extract spike times for this trial
        % Must be a numeric vector
        st = round(spike_times{i});  % round to nearest millisecond

        % Keep only spikes within the window
        st = st(st >= tWin(1) & st <= tWin(2));

        % Create binary spike train (1 ms resolution)
        x = zeros(size(tvec));               % initialize vector of zeros
        x(ismember(tvec, st)) = 1;           % mark timepoints containing spikes

        % Convolve with Gaussian kernel → spike density function
        sdf{i} = conv(x, k, 'same') * Fs;    % convert to spikes/sec
    end
end