function [csType,ipt,saccProps] = detect_catchupSaccade(eyeTraces,stimOnset,pursuitOnset,motionDir,preint,postint,accThresh)
% OBJECTIVE:
% determine if and when a catch-up saccade occurs 
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

if isequal(class(eyeTraces),'struct')
    rPos = eyeTraces.RHPos;
    tPos = eyeTraces.THPos;
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

if (sum(abs(aBase)>accThresh)) > 0 %+sum(abs(vBase)>velThresh)) > 0
    [csAcc,catchup] = max(aBase);
    %[cs_maxVel,idx_maxVel] = max(vBase);

    % saccade amplitude
    rng = rVel; 
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

    % Forward or backward?
    rPos_vec =   rPos(ipt(1)+3:ipt(3)-7);
    rTP_vec = tPos(ipt(1)+3:ipt(3)-7);

    [x1,y1] = pol2cart(rTP_vec(1),rPos_vec(1));
    [x2,y2] = pol2cart(rTP_vec(end),rPos_vec(end));
    x2 = x2 - x1; y2 = y2 - y1;
    [t2,~] = cart2pol(x2,y2);

    cs_angle = wrapTo360(rad2deg(t2));
    cs_diff = abs(wrapTo180(motionDir-cs_angle));
    if cs_diff <= 90
        csType = 2;
    else
        csType = 3;
    end

    % Before or after pursuit onset?
%     blah = 0;
%     if any(ipt<stimOnset) || ipt(1) < stimOnset+50
%         csEpoch = 'fixation';
%     elseif ipt(1) >= stimOnset+50 && ipt(1) < 100
%         csEpoch = 'local';
%     elseif ipt(1) >= pursuitOnset && ipt(1) < pursuitOnset+10
%         csEpoch = 'open loop';
%     elseif ipt()
else
    [ipt,saccProps] = deal(nan(1,3)); csType = 1;
end

end

