function [eyeAcc] = calcDerivative_eyeAcceleration(eyeVel,samplingRate)
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

if isequal(class(eyeVel),'double') 
    % Determine the dimensions of the input data
    [numRows, numCols] = size(eyeVel);

    % Calculate HE and VE acceleration using gradient function
    if numRows > numCols % rows
        x = (0:(size(eyeVel,1) - 1)) / samplingRate;
        HEacc = gradient(eyeVel(:,1)', x)';
        VEacc = gradient(eyeVel(:,2)', x)';
    else % cols
        x = (0:(size(eyeVel,2) - 1)) / samplingRate;
        HEacc = gradient(eyeVel(1,:), x);
        VEacc = gradient(eyeVel(2,:), x);
    end
elseif isequal(class(eyeVel),'cell')
    % Determine the dimensions of the input data
    [~, numCols] = size(eyeVel{1});

    % Calculate HE and VE acceleration using gradient function
    if numCols==2 % rows
        x = cellfun(@(q) (0:(size(q,1) - 1)) / samplingRate, eyeVel, 'uni', 0);
        HEacc = cellfun(@(q,r) gradient(r(:,1)', q)', x, eyeVel, 'uni', 0);
        VEacc = cellfun(@(q,r) gradient(r(:,2)', q)', x, eyeVel, 'uni', 0);

        eyeAcc = cellfun(@(h,v) [h, v], HEacc, VEacc, 'uni', 0);
    else % cols
        x = cellfun(@(q) (0:(size(q,2) - 1)) / samplingRate, eyeVel, 'uni', 0);
        HEacc = cellfun(@(q,r) gradient(r(1,:), q), x, eyeVel, 'uni', 0);
        VEacc = cellfun(@(q,r) gradient(r(2,:), q), x, eyeVel, 'uni', 0);

        eyeAcc = cellfun(@(h,v) [h; v], HEacc, VEacc, 'uni', 0);
    end
end

end

