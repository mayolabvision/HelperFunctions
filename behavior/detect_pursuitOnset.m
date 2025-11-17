function [pursuit_onset, pursuit_latency] = detect_pursuitOnset(eyeVel, stimOnset, targSpeed, varargin)
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
addRequired(p, 'eyeVel', @isnumeric);
addRequired(p, 'stimOnset', @isnumeric);
addRequired(p, 'targSpeed', @isnumeric);
addParameter(p, 'PLOT_TRACES', false, @islogical); % Flag to plot traces

% Parse the inputs
parse(p, eyeVel, stimOnset, targSpeed, varargin{:});

% Assign parsed values to variables
eyeVel = p.Results.eyeVel;
stimOnset = p.Results.stimOnset;
targSpeed = p.Results.targSpeed;
PLOT_TRACES = p.Results.PLOT_TRACES;

fixWin = [-50,50]; % aligned to target motion onset

if isnan(stimOnset)
    pursuit_onset = NaN;
    pursuit_latency = NaN;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert eye velocity to polar coordinates
[~, radVel] = cart2pol(eyeVel(1,:), eyeVel(2,:));

%%%%%%%%%%%%%% 1. DETECT PURSUIT ONSET %%%%%%%%%%%%%%
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
for i = stimOnset + 50:length(radVel) - 10
    if radVel(i) > velocity_threshold && all(radVel(i:i + 9) > threshold_velocity)
        pursuit_onset = i;
        pursuit_latency = pursuit_onset - stimOnset;
        break;
    end
end

% %%%%%%%%%%%%%% 2. DETECT MICROSACCADE  %%%%%%%%%%%%%%
% % Check for microsaccades within the pre-stimulus period
% if (sum(abs(eyeAcc(stimOnset-100:stimOnset+csPre)) > msAcc_thresh) + sum(abs(radVel(stimOnset-100:stimOnset+csPre)) > msVel_thresh)) > 0
%     msOffset = find(abs(eyeAcc(stimOnset-100:stimOnset+csPre)) > msAcc_thresh | abs(radVel(stimOnset-100:stimOnset+csPre)) > msVel_thresh == 1, 1, 'last') - 100;
% else
%     msOffset = NaN;
% end
% 
% %%%%%%%%%%%%%% 3. DETECT CATCH-UP SACCADE %%%%%%%%%%%%%%
% % Check for catch-up saccades within the defined window
% if (sum(abs(eyeAcc(stimOnset + csPre:stimOnset + csPost)) > csAcc_thresh) + sum(abs(radVel(stimOnset + csPre:stimOnset + csPost)) > csVel_thresh) + sum(abs(radVel(stimOnset:stimOnset + crossingTime + 5)) > 15)) > 0 || pursuit_latency > (csPost + crossingTime / 2)
%     % Define binary signal (1 if above threshold, 0 if below)
%     binarySignal = ((abs(eyeAcc(stimOnset + csPre : stimOnset + csPost + 500)) > csAcc_thresh) + ...
%                     (abs(radVel(stimOnset + csPre : stimOnset + csPost + 500)) > csVel_thresh)) > 0;
% 
%     % Find where binarySignal is 1
%     oneIdx = find(binarySignal == 1);
% 
%     % Check for the first occurrence of at least 5 consecutive ones
%     csOnset = NaN;  % Default in case no valid onset is found
%     for i = 1:length(oneIdx) - 4
%         if all(binarySignal(oneIdx(i):oneIdx(i) + 4) == 1)  % Check next 5 frames
%             csOnset = oneIdx(i) + csPre + stimOnset;  % Adjust to original indexing
%             break;
%         end
%     end
% 
%     % Compute catch-up saccade details if onset is found
%     if ~isnan(csOnset)
%         [csVelocity, csPeak] = max(radVel(csOnset:csOnset + 200));
%         csPeak = csPeak + csOnset;
% 
%         % Define binary condition (1 if above threshold, 0 if below)
%         binarySignal = ((abs(eyeAcc(csPeak:csPeak + 200)) > csAcc_thresh) + ...
%                         (abs(radVel(csPeak:csPeak + 200)) > csVel_thresh)) > 0;
% 
%         % Find where binarySignal is 0
%         zeroIdx = find(binarySignal == 0);
% 
%         % Check for the first occurrence of at least 5 consecutive zeros
%         csOffset = NaN;  % Default in case no valid offset is found
%         for i = 1:length(zeroIdx) - 4
%             if all(binarySignal(zeroIdx(i):zeroIdx(i) + 4) == 0)  % Check next 5 frames
%                 csOffset = zeroIdx(i) + csPeak;  % Adjust to original indexing
%                 break;
%             end
%         end
% 
%         % Calculate the catch-up saccade angle
%         dx = eyePos(1, csPeak) - eyePos(1, csOnset);
%         dy = eyePos(2, csPeak) - eyePos(2, csOnset);
%         csAngle = mod(atan2d(dy, dx), 360);
% 
%         % Determine if the catch-up saccade is forward or backward
%         angle_diff = min(mod(csAngle - targAngle, 360), mod(targAngle - csAngle, 360));
%         if angle_diff < 90
%             csType = 'forward';
%         else
%             csType = 'backward';
%         end
%         csFlag = 1;
%     else
%         % If no catch-up saccade detected, set to default 'pure' type
%         csFlag = 0;
%         csOnset = NaN; csVelocity = NaN; csPeak = NaN; csOffset = NaN; csAngle = NaN; angle_diff = NaN; csType = 'pure';
%     end
% else
%     % If no catch-up saccade detected, set to default 'pure' type
%     csFlag = 0;
%     csOnset = NaN; csVelocity = NaN; csPeak = NaN; csOffset = NaN; csAngle = NaN; angle_diff = NaN; csType = 'pure';
% end

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

