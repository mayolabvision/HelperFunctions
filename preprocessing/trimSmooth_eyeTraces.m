function eyes_new = trimSmooth_eyeTraces(eye,stimOnset,preint,postint,win)

if length(eye.HEPos) > (stimOnset+postint-1)
    HEPos = smoothdata(eye.HEPos(stimOnset-preint:stimOnset+postint-1),'gaussian',win).';
    VEPos = smoothdata(eye.VEPos(stimOnset-preint:stimOnset+postint-1),'gaussian',win).';
    HEVel = smoothdata(eye.HEVel(stimOnset-preint:stimOnset+postint-1),'gaussian',win).';
    VEVel = smoothdata(eye.VEVel(stimOnset-preint:stimOnset+postint-1),'gaussian',win).';
    
    x = (1:length(HEVel));
    HEAcc = (gradient(HEVel(:)) ./ gradient(x(:)./1000));
    VEAcc = (gradient(VEVel(:)) ./ gradient(x(:)./1000));
    
    eyes_new = {[HEPos,VEPos]; [HEVel,VEVel]; [HEAcc,VEAcc]};
else
    eyes_new = {NaN};
end

% f = figure;
% 
% subplot(2,2,1)
% x = 1:length(HEPos);
% plot(x,eye.HEPos(stimOnset-preint:stimOnset+postint-1),'k-','linewidth',2)
% hold on
% plot(x,HEPos,'m-','linewidth',2)
% 
% subplot(2,2,2)
% plot(x,eye.VEPos(stimOnset-preint:stimOnset+postint-1),'k-','linewidth',2)
% hold on
% plot(x,VEPos,'m-','linewidth',2)
% 
% subplot(2,2,3)
% plot(x,eye.HEVel(stimOnset-preint:stimOnset+postint-1),'k-','linewidth',2)
% hold on
% plot(x,HEVel,'m-','linewidth',2)
% 
% subplot(2,2,4)
% plot(x,eye.VEVel(stimOnset-preint:stimOnset+postint-1),'k-','linewidth',2)
% hold on
% plot(x,VEVel,'m-','linewidth',2)

end