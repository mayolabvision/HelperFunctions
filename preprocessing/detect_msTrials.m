function msFlag = detect_msTrials(eye,stimOnset,preint,postint,accThresh,velThresh)

% smooth velocities
hVel = smoothdata(eye{3},'gaussian',20); 
vVel = smoothdata(eye{4},'gaussian',20);
[~,rVel] = cart2pol(hVel,vVel); 

x = (1:length(rVel));
rAcc = (gradient(rVel(:)) ./ gradient(x(:)./1000));

% Detect saccades occurring in the window [-51 100], too early in the trial
if sum(abs(rAcc(stimOnset-preint:stimOnset+postint))>accThresh | abs(rVel(stimOnset-preint:stimOnset+postint))'>velThresh)
    msFlag = 1;
else
    msFlag = 0;
end

end

