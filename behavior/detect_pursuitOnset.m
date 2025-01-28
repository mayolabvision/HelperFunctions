function [pursuit_onset,rxnTime,msFlag,csFlag] = detect_pursuitOnset(eyePos,eyeVel,stimOnset,jumpDuration,targSpeed,targDir,jumpType,plt)
% OBJECTIVE:
% determine time of pursuit onset 
%
% INPUTS:
% eyeVel = radial eye velocity either in 1D array of within struct, with fieldname 'RHVel'
% stimOnset = time of target motion onset
% minRT = minimum reaction time (RT) 
% maxRT = maximum reaction time (RT)
%
% OUTPUTS:
% pursuitOnset = time in ms of pursuit onset 
% rxnTime = pursuitOnset - stimOnset 

csPre = 60;
csPost = 210;

% Design a 2nd-order low-pass Butterworth filter
[b, a] = butter(2, 20 / (1000 / 2), 'low');

% Apply the filter to the data
filtered_eyePos = filtfilt(b, a, eyePos);
eyeVel2 = diff(filtered_eyePos) * 1000; % overly smoothed, used to find pursuit onset

% Compute the acceleration using the 10ms before and after difference
window_size = round(10e-3 * 1000);  % Convert 10ms to number of samples

% Initialize the acceleration array
eyeAcc = NaN(1, length(eyeVel2));  % NaN for boundary points where no acceleration can be computed

% Compute acceleration by dividing the velocity difference over 20ms
for i = (window_size + 1):(length(eyeVel2) - window_size)
    eyeAcc(i) = (eyeVel2(i + window_size) - eyeVel2(i - window_size)) / (20e-3);  % 20ms difference
end

% Calculate the standard deviation of the velocity during the fixation period
fixation_velocity = eyeVel2(stimOnset - 50:stimOnset + 50);
fixation_velocity = fixation_velocity(abs(eyeAcc(stimOnset-50:stimOnset+50))<200);
fixation_std = std(fixation_velocity);

% Define the threshold: two times the standard deviation of fixation velocity
velocity_threshold = 2 * fixation_std;

% Search for the pursuit onset: first time velocity > 2 * std(fixation velocity)
% and stays above 30% of target velocity for 5 consecutive frames (10 ms)
threshold_velocity = 0.3 * targSpeed;
pursuit_onset = NaN;  % Initialize pursuit onset index as NaN

for i = stimOnset+50:length(eyeVel2)-10
    if eyeVel2(i) > velocity_threshold && all(eyeVel2(i:i+9) > threshold_velocity)
        pursuit_onset = i;
        break;  % Exit the loop as soon as the condition is satisfied
    end
end

if (sum(abs(eyeAcc(stimOnset-25:stimOnset+csPre))>1000) + sum(abs(eyeVel(stimOnset-25:stimOnset+csPre))>10)) > 0
    msFlag = 1;
    rxnTime = NaN;
    csFlag = NaN;
else
    msFlag = 0;
    rxnTime = pursuit_onset - stimOnset;

    if (sum(abs(eyeAcc(stimOnset+csPre:stimOnset+csPost))>1000) + sum(abs(eyeVel(stimOnset+csPre:stimOnset+csPost))>targSpeed*2)) > 0
        csFlag = 1;
    else
        csFlag = 0;
        if rxnTime > csPost
            csFlag = 1;
        end
    end
end

if plt==1
    fig = figure;
    fig.Position = [100 100 900 900];
    tl = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'loose');

    nexttile

    x = -100:500;
    plot(x,eyePos(stimOnset-100:stimOnset+500),'k-','linewidth',2);


    
    l(1) = xline(0,'k--','linewidth',2);
    hold on
    l(2) = xline(jumpDuration,'k--','linewidth',2);
    if ~isnan(rxnTime)
        l(3) = xline(rxnTime,'b--','linewidth',2);
    else
        l(3) = xline(0, 'b--', 'linewidth', 2);
    end

    %plot(x,eyePos(stimOnset-100:stimOnset+500),'k-','linewidth',2);

    yLimits = ylim; % Get current y-axis limits
    l(4) = fill([csPre csPost csPost csPre], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

    title(sprintf('RT = %d ms, csFlag = %d',rxnTime, csFlag))
    subtitle(sprintf('Speed = %d deg/s, Direction = %d deg, Jump = %d',targSpeed,targDir,jumpType))
    ylabel('position')
    
    legend(l,"step ramp onset", "step ramp offset", "pursuit onset", "detection window",'location','best')
    prettyFig;

    nexttile
    xline(0,'k--','linewidth',2)
    hold on
    xline(jumpDuration,'k--','linewidth',2)
    if ~isnan(rxnTime)
        xline(rxnTime,'b--','linewidth',2)
    end
    yline(targSpeed,'k-','linewidth',2)
    plot(x,eyeVel2(stimOnset-100:stimOnset+500),'k-','linewidth',2)
    ylabel('velocity')
    prettyFig;

    nexttile
    xline(0,'k--','linewidth',2)
    hold on
    xline(jumpDuration,'k--','linewidth',2)
    if ~isnan(rxnTime)
        xline(rxnTime,'b--','linewidth',2)
    end
    plot(x,eyeAcc(stimOnset-100:stimOnset+500),'k-','linewidth',2)
    xlabel('time aligned to target motion onset (ms)')
    ylabel('acceleration')
    prettyFig;
end

end

