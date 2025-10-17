function [eye_vel, eye_acc] = calcDerivative_eyeTraces(EYE, Fs)
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
        Fs = 1000;
    end
    
    % Determine the dimensions of the input data
    if size(EYE,2)>size(EYE,1)
        eye = EYE;
    else
        eye = EYE';
    end
    
    % VELOCITY
    eye_vel = nan(size(eye));
    x = (0:(size(eye,2) - 1)) / Fs;
    for r = 1:size(eye,1)
        eye_vel(r, :) = gradient(eye(r, :), x);
    end

    % ACCELERATION
    dt = 1 / Fs;                     % time per sample (s)
    offset = round(0.010 * Fs);      % number of samples in 10 ms
    
    eye_acc = nan(size(eye_vel));
    for r = 1:size(eye_vel,1)
        % Use central difference: (v[t+10ms] - v[t-10ms]) / 20ms
        eye_acc(r, :) = (circshift(eye_vel(r, :), -offset) - circshift(eye_vel(r, :), offset)) / (2 * offset * dt);
    end

    if size(EYE,2)<size(EYE,1)
        eye_vel = eye_vel';
        eye_acc = eye_acc';
    end

end

