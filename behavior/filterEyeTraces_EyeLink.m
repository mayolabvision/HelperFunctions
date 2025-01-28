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
    
    defaultFs        =  1000; % 
    defaultFc        =  84;
    defaultPlotTrial =  false;
    
    % Create an input parser
    p = inputParser;
    addRequired(p, 'EYE', @isnumeric); % PATH must be a string
    addParameter(p, 'SAMPLING_FREQUENCY', defaultFs, @(x) isnumeric(x));
    addParameter(p, 'CUTOFF_FREQUENCY', defaultFc, @(x) isnumeric(x));
    addParameter(p, 'PLOT_TRIAL', defaultPlotTrial, @islogical); % CORRECT_ONLY must be logical

    % Parse the inputs
    parse(p, EYE, varargin{:});

    % Assign parsed values to variables
    EYE = p.Results.EYE;
    Fs = p.Results.SAMPLING_FREQUENCY;
    Fc = p.Results.CUTOFF_FREQUENCY;
    plot_trial = p.Results.PLOT_TRIAL;
    Ntaps = 80;        % Number of taps (filter order)
    
    % Normalize the cutoff frequency (since fir1 expects normalized frequencies)
    Wn = Fc / (Fs / 2);  % Normalized frequency
    
    % Design FIR filter using Hamming window
    fir_coefficients = fir1(Ntaps - 1, Wn, hamming(Ntaps));  % Use Hamming window for better performance
    
    eye_new = nan(size(EYE));
    if max(size(EYE))>237
        if size(EYE,1) < size(EYE,2)  % EYE is [2, time_points], filter along the time dimension (columns)
            % if size()
            for v=1:size(EYE,1)
                eye_new(v,:) = filtfilt(fir_coefficients, 1, EYE(v,:)); 
            end
    
        else  % EYE is [time_points, 2], filter along the time dimension (rows)
            for v=1:size(EYE,2)
                eye_new(:,v) = filtfilt(fir_coefficients, 1, EYE(:,v));  
            end
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
    
end

