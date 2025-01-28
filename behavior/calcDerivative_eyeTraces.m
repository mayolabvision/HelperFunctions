function [eye_prime] = calcDerivative_eyeTraces(EYE,samplingRate)
% OBJECTIVE:
% Use gradient function to approximate derivative (e.g. eye acceleration)
% Recommended that the eye velocity is smoothed with the sglolay fxn
%
% INPUTS:
% eyeVel = HE/VE eye velocity traces in either a 2D double array or cells of 2D double arrays
%         can be a (2xN), (Nx2) or cells where eyeVel{1} is either (2xN) or (Nx2)
% samplingRate = sampling rate of data, to determine a time vector for appropriate units
%
% OUTPUTS:
% eyeAcc = smoothed eye traces, with same output dimensions as input
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2
        samplingRate = 1000;
    end
    
    % Determine the dimensions of the input data
    [numRows, numCols] = size(EYE);
    
    eye_prime = zeros(size(EYE));
    % Calculate HE and VE acceleration using gradient function
    if numRows > numCols % rows
        x = (0:(size(EYE,1) - 1)) / samplingRate;
        for v=1:size(EYE,2)
            eye_prime(:,v) = gradient(EYE(:,v)', x)';
        end
    else % cols
        x = (0:(size(EYE,2) - 1)) / samplingRate;
        for v=1:size(EYE,1)
            eye_prime(v,:) = gradient(EYE(v,:), x);
        end
    end

end

