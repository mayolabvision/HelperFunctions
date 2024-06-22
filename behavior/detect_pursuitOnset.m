function [pursuitOnset,rxnTime] = detect_pursuitOnset(eyeVelocity,stimOnset,minRT,maxRT,plt)
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

if nargin < 3
    minRT = 50; maxRT = 300;
elseif nargin < 5
    plt = 0;
end

if isequal(class(eyeVelocity),'struct')
    eyeVelocity = eyeVelocity.RHVel;
end

vBase = eyeVelocity(stimOnset-75:stimOnset+25);
if mean(vBase)+(std(vBase)*4) > 5 % if there is a microsaccade, look for a "better" window
    test = eyeVelocity(stimOnset-200:stimOnset+50);
    windows = buffer(test,100,100-10,'nodelay');
    [~, minMeanIndex] = min(mean(windows));
    minMeanStartIndex = (minMeanIndex - 1) * 10 + 1;
    vInd = (stimOnset-200)+minMeanStartIndex;

    vBase = eyeVelocity(vInd:vInd+100); % best window, in terms of vel
end

baseVel = mean(vBase); % baseline eye velocity 
baseVelstd = std(vBase); % STD of baseline eye velocity
stdsBase = (eyeVelocity - baseVel)./baseVelstd; % for each time point, calculate stddev from baseline velocity
stdsBase(1:stimOnset+minRT) = NaN;
stdsBase(stimOnset+maxRT:end) = [];

pursuitOnset = find((stdsBase > baseVel+(baseVelstd*4))==0,1,'last') - 1; % Find last instance eye velocity doesn't exceeds 4 standard deviations
rxnTime = pursuitOnset - stimOnset;

if isempty(rxnTime) || rxnTime >= maxRT-3
    pursuitOnset = find((stdsBase > baseVel+(baseVelstd*3))==0,1,'last') - 1; 
    rxnTime = pursuitOnset - stimOnset;
    if isempty(pursuitOnset) || rxnTime >= maxRT-2
        pursuitOnset = find((stdsBase > baseVel+(baseVelstd*2))==0,1,'last') - 1;  
        if isempty(pursuitOnset) || rxnTime >= maxRT-1
            pursuitOnset = find((stdsBase > baseVel+(baseVelstd*1))==0,1,'last') - 1; 
            plt = 1;
        end
    end
    if ~isempty(pursuitOnset)
        rxnTime = pursuitOnset - stimOnset;
    else
        rxnTime = NaN; pursuitOnset = NaN; 
        return
    end
elseif rxnTime <= minRT+1
    pursuitOnset = find((stdsBase > baseVel+(baseVelstd*5))==0,1,'last') - 1;
    rxnTime = pursuitOnset - stimOnset;
    if rxnTime <= minRT
        pursuitOnset = find((stdsBase > baseVel+(baseVelstd*6))==0,1,'last') - 1;
        rxnTime = pursuitOnset - stimOnset;
        if rxnTime <= minRT
            pursuitOnset = find((stdsBase > baseVel+(baseVelstd*7))==0,1,'last') - 1;
            rxnTime = pursuitOnset - stimOnset;
        end
    end
end

if rxnTime >= maxRT-1 || rxnTime < minRT
    rxnTime = NaN; pursuitOnset = NaN;
    %plt = 1;
end

if plt ==1
    figure;
    x = 1:1500+1;
    plot(x,eyeVelocity(stimOnset-500:stimOnset+1000),'k-')
    hold on
    xline(rxnTime,'b--')
    title(rxnTime)
end



end

