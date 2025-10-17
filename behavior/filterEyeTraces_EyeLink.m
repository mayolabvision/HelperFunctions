function [eye_new] = filterEyeTraces_EyeLink(EYE,varargin)
% OBJECTIVE:
% The Savitzky-Golay (SG) filter, used as digital filter for smoothing noisy data
% Used here for smoothing eye traces
%
% INPUTS:
% EYE = eye traces (1D or 2D) in either a double array or cells of double arrays
% order = polynomial order, choice should balance noise reduction with preservation of important signal features
%         higher order - can fit more intricate shapes, can overfit noise
%         lower order  - simpler fit, can be more robust to noise but less responsive to rapid changes in data
%    (typical values are between 2-5)
% frameLength = time scale over which the filter operates, choice should consider the characteristic time scales of the signal
%         longer frame  - captures long-term trends in data, good for slow-changing signals
%         shorter frame - more responsive to short-term variations, smooths out long-term trends
%    (must be odd, typical is around 9-19)
% plt = whether you want to plot raw & smoothed eye traces
%         0 = don't plot anything
%         1 = plot raw & smoothed traces
%
% OUTPUTS:
% eye_new = smoothed eye traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % Create an input parser
    p = inputParser;
    addRequired(p, 'EYE', @isnumeric); % PATH must be a string
    addParameter(p, 'FILT_TYPE', 'butter', @ischar);
    addParameter(p, 'SAMPLING_FREQUENCY', 1000, @isnumeric);
    addParameter(p, 'CUTOFF_FREQUENCY', 30, @isnumeric);
    addParameter(p, 'ORDER', 4, @isnumeric);
    addParameter(p, 'FRAME_LENGTH', 19, @isnumeric);
    addParameter(p, 'N_TAPS', 80, @isnumeric);
    addParameter(p, 'PLOT_TRIAL', false, @islogical); 

    % Parse the inputs
    parse(p, EYE, varargin{:});

    % Assign parsed values to variables
    EYE = p.Results.EYE;
    filt_type = p.Results.FILT_TYPE;
    Fs = p.Results.SAMPLING_FREQUENCY;
    Fc = p.Results.CUTOFF_FREQUENCY;
    order = p.Results.ORDER;
    frame_length = p.Results.FRAME_LENGTH;
    Ntaps = p.Results.N_TAPS;
    plot_trial = p.Results.PLOT_TRIAL;

    if size(EYE,2)>size(EYE,1)
        eye = EYE;
    else
        eye = EYE';
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
        
        plot(x, EYE, 'k', 'LineWidth', 1, 'DisplayName', 'raw');
        hold on;
        plot(x, eye_new, 'b-',  'LineWidth', 1, 'DisplayName', 'SG filtered');
        xlabel('time aligned to trial onset (ms)');
        ylabel('eye position (deg)');
        legend('location','best');
        prettyFig;
    end

    if size(EYE,2)<size(EYE,1)
        eye_new = eye_new';
    end
    
end

