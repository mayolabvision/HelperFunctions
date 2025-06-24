%% GRC FIGURES

data_path = '/Users/kendranoneman/Data/dualhemi_unleashed';

%% individual session
sess = 'Ya_250515_s395_g0';

load(fullfile(data_path,[sess '_unleashed.mat']), 'S')

%% figure out pursuit onset for Walter

Tpurs = S.pursuit_task_0001.tbl;
Tpurs = Tpurs(Tpurs.result=='CORRECT',:);

trl = 1;
[pursuit_onset, rxnTime, msOffset, csOnset, csVelocity, csPeak, csOffset, csAngle, csType] = detect_pursuitOnset(Tpurs.eyePos{trl}, Tpurs.eyeVel{trl},  Tpurs.PURSUIT_TARG_ON(trl), Tpurs.CROSSING_TIME(trl), Tpurs.pursuitSpeed(trl),  Tpurs.angle(trl));