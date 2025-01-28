function [eye_new] = trimEyeTraces(eye,stimOnset,preint,postint)
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

%preint = stimOnset-1;
%postint = length(eye.HEPos)-stimOnset-1;

if isequal(class(eye),'struct')
    eye_new = structfun(@(q) q(stimOnset-preint:stimOnset+postint), eye, 'uni',0);
elseif isequal(class(eye),'double')
    blah = 0; 
elseif isequal(class(eye),'cell')
    blah = 0;
end

end

