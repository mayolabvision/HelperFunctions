function [pursuit_onset, rxnTime, msOffset, csOnset, csVelocity, csPeak, csOffset, csAngle, csType] = detect_pursuitOnset(eyePos, eyeVel, stimOnset, crossingTime, targSpeed, targAngle, varargin)
% detect_pursuitOnset
% 
% This function detects the onset of pursuit, reaction time, microsaccades,
% and catch-up saccades during a target motion experiment. It also computes 
% the corresponding angular and velocity-related measures.
%
% INPUTS:
%   eyePos       -  2xN matrix of eye position (in Cartesian coordinates)
%   eyeVel       -  2xN matrix of eye velocity (in Cartesian coordinates)
%   stimOnset    -  Time (in ms) when the stimulus (target) onset occurs
%   crossingTime -  Time (in ms) when the target crosses a predefined position
%   targSpeed    -  Speed (in deg/s) of the target motion
%   targAngle    -  Angle (in degrees) of the target motion (relative to some reference)
%
%%%% Optional parameters: %%%
%   'CS_PREINT'   - Pre-catch-up saccade interval (default: 50 ms)
%   'CS_POSTINT'  - Post-catch-up saccade interval (default: 100 ms)
%   'PLOT_TRACES' - Logical flag to plot traces (default: false)
%
% OUTPUTS:
%   pursuit_onset   - Time (in ms) when pursuit onset is detected
%   rxnTime         - Reaction time (in ms) = pursuit_onset - stimOnset
%   msOffset        - Time (in ms) of the last microsaccade before pursuit onset
%   csOnset         - Time (in ms) of the catch-up saccade onset
%   csVelocity      - Peak velocity (in deg/s) of the catch-up saccade
%   csPeak          - Time (in ms) of the peak velocity of the catch-up saccade
%   csOffset        - Time (in ms) of the catch-up saccade offset
%   csAngle         - Angle (in degrees) of the catch-up saccade relative to the target
%   csType          - Type of catch-up saccade ('forward', 'backward', or 'pure')

% Set up input parser
p = inputParser;
addRequired(p, 'eyePos', @isnumeric);
addRequired(p, 'eyeVel', @isnumeric);
addRequired(p, 'stimOnset', @isnumeric);
addRequired(p, 'crossingTime', @isnumeric);
addRequired(p, 'targSpeed', @isnumeric);
addParameter(p, 'CS_PREINT', 50, @isnumeric);  % Pre-catch-up interval (ms)
addParameter(p, 'CS_POSTINT', 100, @isnumeric); % Post-catch-up interval (ms)
addParameter(p, 'PLOT_TRACES', false, @islogical); % Flag to plot traces

% Parse the inputs
parse(p, eyePos, eyeVel, stimOnset, crossingTime, targSpeed, varargin{:});

% Assign parsed values to variables
eyePos = p.Results.eyePos;
eyeVel = p.Results.eyeVel;
stimOnset = p.Results.stimOnset;
crossingTime = p.Results.crossingTime;
targSpeed = p.Results.targSpeed;
CS_PREINT = p.Results.CS_PREINT;
CS_POSTINT = p.Results.CS_POSTINT;
PLOT_TRACES = p.Results.PLOT_TRACES;

%%%%%%%%%%%%%%% BOUNDARIES AND DEFINITIONS %%%%%%%%%%%%%%%%%
csPre = crossingTime - CS_PREINT;  % Pre-catch-up saccade window
csPost = crossingTime + CS_POSTINT;  % Post-catch-up saccade window

msAcc_thresh = 1000; % Microsaccade acceleration threshold (deg/s^2)
msVel_thresh = 10;   % Microsaccade velocity threshold (deg/s)

csAcc_thresh = 1000;  % Catch-up saccade acceleration threshold (deg/s^2)
csVel_thresh = targSpeed * 2;  % Catch-up saccade velocity threshold (deg/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert eye position and velocity to polar coordinates
[~, radPos] = cart2pol(eyePos(1,:), eyePos(2,:));
[~, radVel] = cart2pol(eyeVel(1,:), eyeVel(2,:));

% Design a 2nd-order low-pass Butterworth filter
[b, a] = butter(2, 20 / (1000 / 2), 'low');  % Cutoff at 20 Hz

% Apply the filter to the eye position data
filtered_radPos = filtfilt(b, a, radPos);

% Compute the velocity (first derivative of position)
radVel2 = diff(filtered_radPos) * 1000;  % Convert to deg/s

% Compute acceleration using the 10ms before and after difference
window_size = round(10e-3 * 1000);  % 10 ms in samples

% Initialize the acceleration array
eyeAcc = NaN(1, length(radVel2));  % NaN for boundary points

% Compute acceleration
for i = (window_size + 1):(length(radVel2) - window_size)
    eyeAcc(i) = (radVel2(i + window_size) - radVel2(i - window_size)) / (20e-3);  % 20ms difference
end

%%%%%%%%%%%%%% 1. DETECT PURSUIT ONSET %%%%%%%%%%%%%%
% Calculate the standard deviation of the velocity during the fixation period
fixation_velocity = radVel2(stimOnset - 50:stimOnset + 50);
fixation_velocity = fixation_velocity(abs(eyeAcc(stimOnset-50:stimOnset+50)) < 200);
fixation_std = std(fixation_velocity);

% Define the threshold: two times the standard deviation of fixation velocity
velocity_threshold = 2 * fixation_std;

% Search for the pursuit onset
threshold_velocity = 0.3 * targSpeed;  % 30% of target speed
pursuit_onset = NaN;  % Initialize pursuit onset index
rxnTime = NaN;

% Loop through the data to find pursuit onset
for i = stimOnset + 50:length(radVel2) - 10
    if radVel2(i) > velocity_threshold && all(radVel2(i:i + 9) > threshold_velocity)
        pursuit_onset = i;
        rxnTime = pursuit_onset - stimOnset;
        break;
    end
end

%%%%%%%%%%%%%% 2. DETECT MICROSACCADE  %%%%%%%%%%%%%%
% Check for microsaccades within the pre-stimulus period
if (sum(abs(eyeAcc(stimOnset-100:stimOnset+csPre)) > msAcc_thresh) + sum(abs(radVel(stimOnset-100:stimOnset+csPre)) > msVel_thresh)) > 0
    msOffset = find(abs(eyeAcc(stimOnset-100:stimOnset+csPre)) > msAcc_thresh | abs(radVel(stimOnset-100:stimOnset+csPre)) > msVel_thresh == 1, 1, 'last') - 100;
else
    msOffset = NaN;
end

%%%%%%%%%%%%%% 3. DETECT CATCH-UP SACCADE %%%%%%%%%%%%%%
% Check for catch-up saccades within the defined window
if (sum(abs(eyeAcc(stimOnset + csPre:stimOnset + csPost)) > csAcc_thresh) + sum(abs(radVel(stimOnset + csPre:stimOnset + csPost)) > csVel_thresh)) > 0 || rxnTime > (csPost + crossingTime / 2)
    % Define binary signal (1 if above threshold, 0 if below)
    binarySignal = ((abs(eyeAcc(stimOnset + csPre : stimOnset + csPost + 500)) > csAcc_thresh) + ...
                    (abs(radVel(stimOnset + csPre : stimOnset + csPost + 500)) > csVel_thresh)) > 0;

    % Find where binarySignal is 1
    oneIdx = find(binarySignal == 1);

    % Check for the first occurrence of at least 5 consecutive ones
    csOnset = NaN;  % Default in case no valid onset is found
    for i = 1:length(oneIdx) - 4
        if all(binarySignal(oneIdx(i):oneIdx(i) + 4) == 1)  % Check next 5 frames
            csOnset = oneIdx(i) + csPre + stimOnset;  % Adjust to original indexing
            break;
        end
    end

    % Compute catch-up saccade details if onset is found
    if ~isnan(csOnset)
        [csVelocity, csPeak] = max(radVel(csOnset:csOnset + 200));
        csPeak = csPeak + csOnset;

        % Define binary condition (1 if above threshold, 0 if below)
        binarySignal = ((abs(eyeAcc(csPeak:csPeak + 200)) > csAcc_thresh) + ...
                        (abs(radVel(csPeak:csPeak + 200)) > csVel_thresh)) > 0;

        % Find where binarySignal is 0
        zeroIdx = find(binarySignal == 0);

        % Check for the first occurrence of at least 5 consecutive zeros
        csOffset = NaN;  % Default in case no valid offset is found
        for i = 1:length(zeroIdx) - 4
            if all(binarySignal(zeroIdx(i):zeroIdx(i) + 4) == 0)  % Check next 5 frames
                csOffset = zeroIdx(i) + csPeak;  % Adjust to original indexing
                break;
            end
        end

        % Calculate the catch-up saccade angle
        dx = eyePos(1, csPeak) - eyePos(1, csOnset);
        dy = eyePos(2, csPeak) - eyePos(2, csOnset);
        csAngle = mod(atan2d(dy, dx), 360);

        % Determine if the catch-up saccade is forward or backward
        angle_diff = min(mod(csAngle - targAngle, 360), mod(targAngle - csAngle, 360));
        if angle_diff < 90
            csType = 'forward';
        else
            csType = 'backward';
        end
    else
        % If no catch-up saccade detected, set to default 'pure' type
        csFlag = 0;
        csOnset = NaN; csVelocity = NaN; csPeak = NaN; csOffset = NaN; csAngle = NaN; angle_diff = NaN; csType = 'pure';
    end
else
    % If no catch-up saccade detected, set to default 'pure' type
    csFlag = 0;
    csOnset = NaN; csVelocity = NaN; csPeak = NaN; csOffset = NaN; csAngle = NaN; angle_diff = NaN; csType = 'pure';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PLOT_TRACES
    fig = figure;
    fig.Position = [100 100 1300 900];
    tl = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'loose');

    x = -25:500;

    nexttile  
    plot(x,radPos(stimOnset-25:stimOnset+500),'k-','linewidth',2);

    l(1) = xline(0,'k--','linewidth',2);
    hold on
    l(2) = xline(crossingTime,'k--','linewidth',2);
    if ~isnan(rxnTime)
        l(3) = xline(rxnTime,'b--','linewidth',2);
    else
        l(3) = xline(0, 'b--', 'linewidth', 2);
    end

    %plot(x,radPos(stimOnset-100:stimOnset+500),'k-','linewidth',2);

    yLimits = ylim; % Get current y-axis limits
    l(4) = fill([csPre csPost csPost csPre], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

    if isequal(csType,'forward')
        xline(csOnset-stimOnset,'g--','linewidth',2)
        xline(csPeak-stimOnset,'g--','linewidth',2)
        xline(csOffset-stimOnset,'g--','linewidth',2)
    elseif isequal(csType,'backward')
        xline(csOnset-stimOnset,'r--','linewidth',2)
        xline(csPeak-stimOnset,'r--','linewidth',2)
        xline(csOffset-stimOnset,'r--','linewidth',2)
    end
    title(sprintf('RT = %d ms, csFlag = %d',rxnTime, csFlag))
    subtitle(sprintf('Speed = %d deg/s, Angle = %d deg',targSpeed,targAngle))
    ylabel('position')
    xlim([-25,500])
    legend(l,"step ramp onset", "step ramp offset", "pursuit onset", "detection window",'location','best')
    prettyFig;

    nexttile
    xline(0,'k--','linewidth',2)
    hold on
    xline(crossingTime,'k--','linewidth',2)
    if ~isnan(rxnTime)
        xline(rxnTime,'b--','linewidth',2)
    end
    if isequal(csType,'forward')
        xline(csOnset-stimOnset,'g--','linewidth',2)
        xline(csPeak-stimOnset,'g--','linewidth',2)
        xline(csOffset-stimOnset,'g--','linewidth',2)
    elseif isequal(csType,'backward')
        xline(csOnset-stimOnset,'r--','linewidth',2)
        xline(csPeak-stimOnset,'r--','linewidth',2)
        xline(csOffset-stimOnset,'r--','linewidth',2)
    end
    yline(targSpeed,'k-','linewidth',2)
    plot(x,radVel2(stimOnset-25:stimOnset+500),'k-','linewidth',2)
    ylabel('velocity')
    xlim([-25,500])
    prettyFig;

    nexttile
    xline(0,'k--','linewidth',2)
    hold on
    xline(crossingTime,'k--','linewidth',2)
    if ~isnan(rxnTime)
        xline(rxnTime,'b--','linewidth',2)
    end
    plot(x,eyeAcc(stimOnset-25:stimOnset+500),'k-','linewidth',2)
    xlabel('time aligned to target motion onset (ms)')
    ylabel('acceleration')
    xlim([-25,500])
    prettyFig;
end

end

