function [eye_new] = filterEyeTraces_EyeLink(eyedata,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE:
%   Apply one of three smoothing/filtering methods to raw EyeLink eye traces
%   (x- and y-position). This function provides a unified interface for:
%       • Butterworth low-pass filtering
%       • Savitzky–Golay polynomial smoothing
%       • FIR low-pass filtering
%   The goal is to remove high-frequency noise while preserving true eye
%   movement structure (saccades, smooth pursuit, fixations).
%
% =========================================================================
%
% INPUTS:
%   eyedata  = raw eye traces in a numeric matrix
%              Size should be [2 x T] or [T x 2] for horizontal & vertical
%              eye traces; T = number of time points 
%
%   FILT_TYPE (string) = type of smoothing filter to apply.
%       Options:
%           'butter' – Butterworth low-pass filter (default)
%           'sg'     – Savitzky–Golay filter
%           'fir'    – Windowed FIR low-pass filter
%
%   SAMPLING_FREQUENCY (Hz) = sampling rate of the eye data 
%       Default: 1000 Hz, for EyeLink eye tracker
%
%   PLOT_TRIAL (boolean) = whether to plot raw vs filtered traces.
%       0 or false – no plotting (default)
%       1 or true  – plot raw + filtered signals for debugging
%
% -------------------------------------------------------------------------
%
% FILTER-SPECIFIC PARAMETERS
%
% **1. BUTTERWORTH LOW-PASS FILTER ('butter')**
%       - A smooth low-pass filter with a maximally flat frequency response.
%       - Good when you want to attenuate high-frequency noise while preserving
%         overall signal shape.
%
%   CUTOFF_FREQUENCY (Hz)
%       Meaning: frequency above which signal components are attenuated.
%       Default: 30 Hz
%
%
% **2. SAVITZKY–GOLAY FILTER ('sg')**
%       - Fits a sliding polynomial to the data via least squares.
%       - Good for smoothing while preserving peaks, derivatives, and saccades.
%       - Does *not* behave like a frequency-based filter; instead, a local
%         polynomial approximation.
%
%   ORDER
%       Meaning: polynomial order, choice should balance noise reduction with 
%                preservation of important signal features.
%       Default: 4
%       Notes: 
%           Higher order - can fit more intricate shapes, can overfit noise
%           Lower order  - simpler fit, can be more robust to noise but less
%
%   FRAME_LENGTH
%       Meaning: size of the moving window (must be odd).
%       Default: 19 samples
%       Typical: 9–19 samples
%       Notes:
%           Longer windows capture slow trends; shorter windows preserve rapid
%           changes like saccades.
%
%
% **3. FIR LOW-PASS FILTER ('fir')**
%       - A linear-phase filter (preserves waveform shape), using Hamming window.
%       - Useful when phase distortion must be minimized (but filtfilt() already
%         enforces zero-phase).
%
%   N_TAPS
%       Meaning: number of filter coefficients (“length” of the FIR filter).
%       Default: 80
%       Notes:
%           Larger N_TAPS = smoother filter with sharper cutoff, but slower.
%           Rule of thumb: choose N ~ sampling_rate / desired_cutoff.
%
% =========================================================================
%
% OUTPUTS:
%   eye_new = filtered eye traces (same shape as input)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% Create an input parser
p = inputParser;
addRequired(p, 'eyedata', @isnumeric); 
addParameter(p, 'FILT_TYPE', 'butter', @ischar); % 'butter', 'sg', or 'fir'
addParameter(p, 'SAMPLING_FREQUENCY', 1000, @isnumeric);
addParameter(p, 'CUTOFF_FREQUENCY', 30, @isnumeric);
addParameter(p, 'ORDER', 4, @isnumeric); % only used for SGL filter
addParameter(p, 'FRAME_LENGTH', 19, @isnumeric); % only used for SGL filter
addParameter(p, 'N_TAPS', 80, @isnumeric); % only used for FIR filter
addParameter(p, 'PLOT_TRIAL', false, @islogical); 

% Parse the inputs
parse(p, eyedata, varargin{:});

% Assign parsed values to variables
eyedata = p.Results.eyedata;
filt_type = p.Results.FILT_TYPE;
Fs = p.Results.SAMPLING_FREQUENCY;
Fc = p.Results.CUTOFF_FREQUENCY;
order = p.Results.ORDER;
frame_length = p.Results.FRAME_LENGTH;
Ntaps = p.Results.N_TAPS;
plot_trial = p.Results.PLOT_TRIAL;

if size(eyedata,2)>size(eyedata,1)
    eye = eyedata;
else
    eye = eyedata';
end

eye_new = nan(size(eye));
if isequal(filt_type, 'butter') % butterworth
    [b, a] = butter(order, Fc / (Fs / 2), 'low');

    for r = 1:size(eye_new,1)
        eye_new(r, :) = filtfilt(b, a, eye(r, :));
    end

elseif isequal(filt_type, 'sg') % savitsky-golay
    for r = 1:size(eye_new,1)
        eye_new(r, :) = sgolayfilt(eye(r, :), order, frame_length); 
    end

elseif isequal(filt_type, 'fir')
    Wn = Fc / (Fs / 2);  % Normalized frequency
    fir_coefficients = fir1(Ntaps - 1, Wn, hamming(Ntaps));

    for r = 1:size(eye_new,1)
        eye_new(r, :) = filtfilt(fir_coefficients, 1, eye(f, :)); 
    end
end

if plot_trial
    figure('Position', [200, 100, 1400, 600]);
    x = (1:length(eye_new));
    
    plot(x, eye, 'k', 'LineWidth', 1, 'DisplayName', 'raw');
    hold on;
    plot(x, eye_new, 'b-',  'LineWidth', 1, 'DisplayName', 'SG filtered');
    xlabel('time aligned to trial onset (ms)');
    ylabel('eye position (deg)');
    legend('location','best');
    prettyFig;
end

if size(eyedata,2)<size(eyedata,1)
    eye_new = eye_new';
end
    
end

