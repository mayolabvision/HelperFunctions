function [csType,ipt,saccProps] = detect_catchupSaccade(eye,pursuitOnset,motionDir,preint,postint,accThresh,velThresh)
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

[tPos,rPos] = cart2pol(eye{1},eye{2});
rVel = eye{6};
rAcc = eye{7};

vBase = rVel; vBase(1:pursuitOnset-preint) = NaN; vBase(pursuitOnset+postint:end) = NaN;
aBase = rAcc; aBase(1:pursuitOnset-preint) = NaN; aBase(pursuitOnset+postint:end) = NaN;

if any(abs(aBase)>accThresh & abs(vBase)>velThresh)
    [csAcc,catchup] = max(aBase);

    % saccade amplitude
    rng = rVel; 
    rng(1:catchup-50) = NaN; rng(catchup+50:end) = NaN;
    [csVel,i] = max(rng); 
    catchup = i;

    cs_start = catchup - (76 - find(abs(rAcc(i-75:i))>accThresh,1,'first'));
    cs_end = catchup + find(abs(rAcc(i:i+75))>accThresh,1,'last');

    if (cs_start-pursuitOnset)<postint && (cs_start-pursuitOnset)>-preint && ~isempty(cs_start) && ~isempty(cs_end)
        r1 = rPos(cs_start); r2 = rPos(cs_end);
        t1 = tPos(cs_start); t2 = tPos(cs_end);
        sacAmp = sqrt(r1^2 + r2^2 - 2*(r1*r2)*cos(t1-t2));

        ipt = [cs_start, catchup, cs_end];
        saccProps = [csAcc, csVel, sacAmp];
    else
        [ipt,saccProps] = deal(nan(1,3)); csType = 0;

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
    
else
    [ipt,saccProps] = deal(nan(1,3)); csType = 1;
end

end

