function [pursuit_onset, pursuit_latency] = detect_pursuitOnset(eyeVelocity, stimOnset, targSpeed, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE:
%   Detect smooth pursuit onset and compute pursuit latency relative to target
%   motion onset in a smooth-pursuit eye movement experiment.
%
%   This function implements a physiologically-motivated pursuit detection 
%   algorithm that:
%       1. Converts 2D (horizontal/vertical) velocity into radial velocity.
%       2. Establishes a noise baseline using pre-motion fixation velocity.
%       3. Computes a velocity threshold based on fixation variability.
%       4. Searches for the earliest time point when eye velocity:
%             • Exceeds baseline noise by 2×STD of fixation velocity, AND
%             • Sustains >30% of the target velocity for at least 10 ms.
% 
%    Ref: Goettker, A., Brenner, E., Gegenfurtner, K.R. et al. 
%         Corrective saccades influence velocity judgments and interception. 
%         Sci Rep 9, 5395 (2019). https://doi.org/10.1038/s41598-019-41857-z
%
% =========================================================================
%
% INPUTS:
%   eyeVelocity  
%       Raw 2D eye velocity traces (deg/s).  
%       Expected shape:  
%           - [2 × T] or [T × 2] for horizontal (X) and vertical (Y) velocity  
%           - T = number of samples  
%       Data should already be differentiated/filtered to produce velocity.
%
%   stimOnset  
%       Time index (in samples, not ms) corresponding to the onset of target  
%       motion. If the sampling rate is 1000 Hz, stimOnset=500 corresponds  
%       to 500 ms.  
%       If stimOnset is NaN, the function returns NaN for both outputs.
%
%   targSpeed  
%       The target motion speed (deg/s).  
%       Used when determining sustained pursuit engagement (>30% of target).
%
% -------------------------------------------------------------------------
%
% OPTIONAL PARAMETERS (Name–Value pairs):
%   PLOT_TRACES (default = false)
%       Logical flag controlling whether pursuit detection diagnostics  
%       are plotted:
%         false → no plotting (default)  
%         true  → plot radial velocity, target velocity, and detected onset  
%
% =========================================================================
% OUTPUTS:
%   pursuit_onset  
%       Sample index of detected pursuit onset (relative to recording start).  
%       If no valid onset is found, returns NaN.
%
%   pursuit_latency  
%       Time difference (in samples) between pursuit onset and stimOnset:  
%           pursuit_latency = pursuit_onset – stimOnset  
%
%       If sampling rate = 1000 Hz, this equals latency in milliseconds.  
%       If onset is not detected, returns NaN.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up input parser
p = inputParser;
addRequired(p, 'eyeVel', @isnumeric);
addRequired(p, 'stimOnset', @isnumeric);
addRequired(p, 'targSpeed', @isnumeric);
addParameter(p, 'PLOT_TRACES', false, @islogical); % Flag to plot traces

% Parse the inputs
parse(p, eyeVelocity, stimOnset, targSpeed, varargin{:});

% Assign parsed values to variables
eyeVelocity = p.Results.eyeVel;
stimOnset = p.Results.stimOnset;
targSpeed = p.Results.targSpeed;
PLOT_TRACES = p.Results.PLOT_TRACES;

fixWin = [-50,50]; % aligned to target motion onset

if isnan(stimOnset)
    pursuit_onset = NaN;
    pursuit_latency = NaN;
    return
end

if size(eyeVelocity,2)<size(eyeVelocity,1)
    eyeVelocity = eyeVelocity';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert eye velocity to polar coordinates
[~, radVel] = cart2pol(eyeVelocity(1,:), eyeVelocity(2,:));

% Calculate the standard deviation of the velocity during the fixation period
fixation_velocity = radVel(stimOnset + fixWin(1) : stimOnset + fixWin(2));
fixation_velocity = fixation_velocity(fixation_velocity<3);
fixation_std = std(fixation_velocity);

% Define the threshold: two times the standard deviation of fixation velocity
velocity_threshold = 2 * fixation_std;

% Search for the pursuit onset
threshold_velocity = 0.3 * targSpeed;  % 30% of target speed
pursuit_onset = NaN;  % Initialize pursuit onset index
pursuit_latency = NaN;

% Loop through the data to find pursuit onset
for i = stimOnset + 50 : length(radVel) - 110   % leave room for 100 future points

    % --- Immediate rise requirement (threshold_velocity for 10 consecutive samples)
    immediate_ok = radVel(i) > velocity_threshold && ...
                   all(radVel(i : i+9) > threshold_velocity);

    if ~immediate_ok
        continue
    end

    % --- Sustained requirement (must remain above 2 for 100 points)
    sustain_window = radVel(i : i + 99);
    sustained_ok = all(sustain_window > 1);

    if ~sustained_ok
        % fails the "long window" condition → keep searching
        continue
    end

    % If both conditions pass → true pursuit onset found
    pursuit_onset = i;
    break
end

if ~isnan(pursuit_onset) 
    pursuit_onset2 = pursuit_onset - (50-find(diff(radVel(pursuit_onset-50:pursuit_onset))<0,1,'last')-1);
    if ~isempty(pursuit_onset2)
        pursuit_onset = pursuit_onset2;
    end

    pursuit_latency = pursuit_onset - stimOnset;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PLOT_TRACES
    fig = figure;
    fig.Position = [100 100 1300 900];
   
    x = -25:500;

    xline(0,'k--','linewidth',2)
    hold on
    if ~isnan(pursuit_latency)
        xline(pursuit_latency,'b--','linewidth',2)
    end
   
    yline(targSpeed,'k-','linewidth',2)
    plot(x,radVel(stimOnset-25:stimOnset+500),'k-','linewidth',2)
    ylabel('radial velocity')
    xlim([-25,500])
    ylim([0,30])
    title(sprintf('RT = %d ms',pursuit_latency))
    prettyFig;

end

end

