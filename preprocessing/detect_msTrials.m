function msFlag = detect_msTrials(eye,stimOnset,preint,postint,accThresh,velThresh)
% OBJECTIVE:
% determine if nhp made microsaccade around stimulus onset in given trial
%
% INPUTS:
% eye = 4x1 cell array w/ (HEPos, VEPos, HEVel, VEVel)
% stimOnset = time of target motion onset or other stimulus onset
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

% smooth velocities
hVel = smoothdata(eye{3},'gaussian',20); 
vVel = smoothdata(eye{4},'gaussian',20);
[~,rVel] = cart2pol(hVel,vVel); 

x = (1:length(rVel));
rAcc = (gradient(rVel(:)) ./ gradient(x(:)./1000));

% Detect saccades occurring in the window [-preint postint], where 0 = stim onset
if sum(abs(rAcc(stimOnset-preint:stimOnset+postint))>accThresh | abs(rVel(stimOnset-preint:stimOnset+postint))'>velThresh)
    msFlag = 1;
else
    msFlag = 0;
end

end

