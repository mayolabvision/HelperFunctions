function [msFlag] = detect_microSaccades(eyeTraces,stimOnset,preint,postint,accThresh,velThresh)
% OBJECTIVE:
% determine if nhp made microsaccade around stimulus onset in given trial
%
% INPUTS:
% eye = either a 2D array or 2 cell arrays (RHVel,RHAcc), or a struct with fieldnames 'RHVel' and 'RHAcc' 
% stimOnset = time of target motion onset or other time to align to
% preint = time before stimOnset you want to include 
% postint = time after stimOnset you want to include
% accThresh = acceleration threshold for microsacc detection (e.g. 750 deg/s^2)
% velThresh = velocity threshold for microsacc detection (e.g. 50 deg/s)
%
% OUTPUTS:
% msFlag = 1 if microsacc detected, 0 if not detected

if nargin < 3
    preint = 50; postint = 50; accThresh = 750; velThresh = 50;
elseif nargin < 5
    accThresh = 750; velThresh = 50;
end

if isequal(class(eyeTraces),'struct')
    rVel = eyeTraces.RHVel;
    rAcc = eyeTraces.RHAcc;
elseif isequal(class(eyeTraces),'cell')
    rVel = eyeTraces{1};
    rAcc = eyeTraces{2};
elseif isequal(class(eyeTraces),'double')
    [~, numCols] = size(eyeTraces);
    if numCols==2
        rVel = eyeTraces(:,1);
        rAcc = eyeTraces(:,2);
    else
        rVel = eyeTraces(1,:);
        rAcc = eyeTraces(2,:);
    end
end

vBase = rVel; vBase(1:stimOnset-preint) = NaN; vBase(stimOnset+postint:end) = NaN;
aBase = rAcc; aBase(1:stimOnset-preint) = NaN; aBase(stimOnset+postint:end) = NaN;

if (sum(aBase>accThresh) + sum(vBase>velThresh)) > 0
    msFlag = 1;
else
    msFlag = 0;
end

end

