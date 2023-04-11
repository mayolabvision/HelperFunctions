function [pursuitOnset,rxnTime] = detect_pursuitOnset(eye,stimOnset,minRT,maxRT)
% OBJECTIVE:
% determine time of pursuit onset 
%
% INPUTS:
% eye = 7x1 cell array w/ (HEPos, VEPos, HEVel, VEVel, THVel, RHVel, RHAcc)
% stimOnset = time of target motion onset
% minRT = minimum reaction time (RT) 
% maxRT = maximum reaction time (RT)
%
% OUTPUTS:
% pursuitOnset = time in ms of pursuit onset 
% rxnTime = pursuitOnset - stimOnset 

if nargin < 3
    minRT = 50; maxRT = 300;
end

rVel = eye{6};

vBase = rVel(stimOnset-50:stimOnset+50);
baseVel = mean(vBase); % baseline eye velocity 
baseVelstd = std(vBase); % STD of baseline eye velocity
stdsBase = (rVel - baseVel)./baseVelstd; % for each time point, calculate stddev from baseline velocity
stdsBase(1:stimOnset+minRT) = NaN;
stdsBase(stimOnset+maxRT:end) = [];

pursuitOnset = find((stdsBase > baseVel+(baseVelstd*4))==0,1,'last') + 1; % Find last instance eye velocity doesn't exceeds 4 standard deviations
rxnTime = pursuitOnset - stimOnset;

if rxnTime>250
    pursuitOnset2 = find((stdsBase > baseVel+(baseVelstd*2))==1,1,'first');
    rxnTime2 = pursuitOnset2-stimOnset;
    if rxnTime2<=250
        pursuitOnset = pursuitOnset2;
        rxnTime = rxnTime2;
    else
        pursuitOnset = NaN; rxnTime = NaN;
        return
    end
end
if isempty(rxnTime) || isempty(pursuitOnset)
    rxnTime = 0; pursuitOnset = 0; 
    return
end

end

