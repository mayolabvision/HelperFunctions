function [eye_new] = filterEyeTraces_sglolay(eye,polynomialOrder,frameLength,plt)
% OBJECTIVE:
% The Savitzky-Golay (SG) filter, used as digital filter for smoothing noisy data
% Used here for smoothing eye traces
%
% INPUTS:
% eye = eye traces (1D or 2D) in either a double array or cells of double arrays
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

if nargin < 2
    polynomialOrder = 2; frameLength = 11; plt = 0;
elseif nargin < 4
    plt = 0;
end

if isequal(class(eye),'double')
    % Determine the dimensions of the input data
    [numRows, numCols] = size(eye);

    % Apply the Savitzky-Golay filter along one dimension 
    if numRows > numCols % rows
        eye_new = sgolayfilt(eye, polynomialOrder, frameLength, [], 1);
    else % cols
        eye_new = sgolayfilt(eye', polynomialOrder, frameLength, [], 1)';
    end
elseif isequal(class(eye),'cell')
    % Determine the dimensions of the input data
    [numRows, numCols] = size(eye{1});

    % Apply the Savitzky-Golay filter along one dimension 
    if numRows > numCols % rows
        eye_new = cellfun(@(q) sgolayfilt(q, polynomialOrder, frameLength, [], 1), eye, 'uni', 0);
    else % cols
        eye_new = cellfun(@(q) sgolayfilt(q', polynomialOrder, frameLength, [], 1)', eye, 'uni', 0);
    end
end

if plt == 1 && isequal(class(eye),'double')
    figure;
    x = (1:length(eye_new));
    
    plot(x, eye, 'k', 'LineWidth', 1, 'DisplayName', 'raw');
    hold on;
    plot(x, eye_new, 'b-',  'LineWidth', 1, 'DisplayName', 'SG filtered');
    plot(x, smoothdata(eye,'gaussian',20), 'r-',  'LineWidth', 1, 'DisplayName', 'smoothed');
    xlabel('X');
    ylabel('Y');
    legend('location','best');
elseif plt == 1 && isequal(class(eye),'cell')
    figure;
    x = (1:length(eye_new{1}));
    
    plot(x, eye{1}, 'k', 'LineWidth', 1, 'DisplayName', 'raw');
    hold on;
    plot(x, eye_new{1}, 'b-',  'LineWidth', 1, 'DisplayName', 'SG filtered');
    xlabel('X');
    ylabel('Y');
    legend('location','best');
end

end

