function sacOnset = detect_saccadeOnset(eyeAcc,stimOnset,coordSystem,minRT,maxRT,toPlot)
% OBJECTIVE:
% determine if and when a catch-up saccade occurs 
%
% INPUTS:
% eye = 7x1 cell array w/ (HEPos, VEPos, HEVel, VEVel, ThAcc, RhAcc, RHAcc)
% stimOnset = time of target motion onset
% coordSystem: 'cart' or 'pol'
% minRT = minimum reaction time (RT) 
% maxRT = maximum reaction time (RT)
%
% OUTPUTS:
% pursuitOnset = time in ms of pursuit onset 
% rxnTime = pursuitOnset - stimOnset 

if nargin < 3
    coordSystem = 'pol';
    minRT = 100;
    maxRT = 400;
    toPlot = 0;
elseif nargin < 4
    minRT = 100;
    maxRT = 400;
    toPlot = 0;
end

if isequal(coordSystem,'pol')
    [~, numCols] = size(eyeAcc);
    if numCols==2
        rAcc = eyeAcc(:,1);
    else
        rAcc = eyeAcc(1,:);
    end
elseif isequal(coordSystem,'cart')
    [~, numCols] = size(eyeAcc);
    if numCols==2
        hAcc = eyeAcc(:,1);
        vAcc = eyeAcc(:,2);
    else
        hAcc = eyeAcc(1,:);
        vAcc = eyeAcc(2,:);
    end
    [~,rAcc] = cart2pol(hAcc,vAcc);
end

baseVel = rAcc(stimOnset-100:stimOnset+100);
rAcc(1:stimOnset+minRT-1) = NaN; rAcc(stimOnset+maxRT+1:end) = NaN;   

vBase = rAcc; vBase(1:stimOnset-preint) = NaN; vBase(stimOnset+postint:end) = NaN;
aBase = rAcc; aBase(1:stimOnset-preint) = NaN; aBase(stimOnset+postint:end) = NaN;

if (sum(abs(aBase)>accThresh)) > 0 %+sum(abs(vBase)>velThresh)) > 0
    [csAcc,catchup] = max(aBase);
    %[cs_maxVel,idx_maxVel] = max(vBase);

    % saccade amplitude
    rng = rAcc; 
    rng(1:catchup-50) = NaN; rng(catchup+50:end) = NaN;
    [csVel,i] = max(rng); 
    catchup = i;

    cs_start = catchup - (76 - find(abs(rAcc(i-75:i))>accThresh,1,'first'));
    cs_end = catchup + find(abs(rAcc(i:i+75))>accThresh,1,'last');

    if ~isempty(cs_start) && ~isempty(cs_end)
        if (cs_start-stimOnset)<postint && (cs_start-stimOnset)>-preint 
            r1 = rPos(cs_start); r2 = rPos(cs_end);
            t1 = tPos(cs_start); t2 = tPos(cs_end);
            sacAmp = sqrt(r1^2 + r2^2 - 2*(r1*r2)*cos(t1-t2));
    
            ipt = [cs_start, catchup, cs_end];
            saccProps = [csAcc, csVel, sacAmp];
        else
            [ipt,saccProps] = deal(nan(1,3)); csType = 4;
    
            return
        end
    else
        [ipt,saccProps] = deal(nan(1,3)); csType = 4;
    
        return
    end

end

