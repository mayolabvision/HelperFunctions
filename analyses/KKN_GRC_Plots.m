%% GRC FIGURES

data_path = '/Users/kendranoneman/Data/dualhemi_unleashed';

%% individual session
sess = 'kendra_scrappy_0142a_g0';

load(fullfile(data_path,[sess '_unleashed.mat']), 'S')

%% adding some additional thresholds
FR_thresh     =  1;    % Hz
seldir_range =  [0 0.99]; 
msOff_thresh  =  -50; % ms
sacLat_range =  [100 300];  % ms

%---- Behavior criteria ----%
% mdir criteria 
Tmdir = S.mdir1.tbl(S.mdir1.tbl.result=='CORRECT',:);
Tmdir.saccadeLatency = Tmdir.SACCADE-cellfun(@(q) q(1), Tmdir.FIX_OFF);
Tmdir = Tmdir(Tmdir.saccadeLatency >= sacLat_range(1) & Tmdir.saccadeLatency < sacLat_range(2),:);

[th,rh] = cellfun(@(q) cart2pol(q(1,:),q(2,:)), Tmdir.eyePos, 'uni', 0);
Tmdir = Tmdir(cellfun(@(r,t,s,a) (max(r(s-10:s+75)) < 30) && (abs(wrapTo180(mean(rad2deg(t(s-10:s+75))) - a)) < 45), rh, th, num2cell(Tmdir.SACCADE), num2cell(Tmdir.angle), 'uni', 1),:);

% purs criteria
Tpurs = S.purs1.tbl(S.purs1.tbl.result=='CORRECT' & S.purs1.tbl.pursType=='pure',:);
Tpurs = Tpurs(isnan(Tpurs.msOffset) | Tpurs.msOffset<msOff_thresh,:);

%---- Cluster criteria ----%
Cleft = S.kilosort(1).clusters; Crght = S.kilosort(2).clusters;

% firing rate
Cleft = Cleft(Cleft.mdir1_Hz > FR_thresh & Cleft.purs1_Hz > FR_thresh,:);
Crght = Crght(Crght.mdir1_Hz > FR_thresh & Crght.purs1_Hz > FR_thresh,:);

% selectivity
Cleft = Cleft((Cleft.vis_sel_dir>seldir_range(1) & Cleft.vis_sel_dir<seldir_range(2)) & (Cleft.sac_sel_dir>seldir_range(1) & Cleft.sac_sel_dir<seldir_range(2)) & (Cleft.pur_sel_dir>seldir_range(1) & Cleft.pur_sel_dir<seldir_range(2)),:);
Crght = Crght((Crght.vis_sel_dir>seldir_range(1) & Crght.vis_sel_dir<seldir_range(2)) & (Crght.sac_sel_dir>seldir_range(1) & Crght.sac_sel_dir<seldir_range(2)) & (Crght.pur_sel_dir>seldir_range(1) & Crght.pur_sel_dir<seldir_range(2)),:);

%% behavior distributions

%f = figure;
%histStyle_KKN(Tmdir.saccadeLatency)

% f = figure;
% G = groupsummary(Tmdir, 'delay', {'mean', 'std'}, 'saccadeLatency');
% errorbar(G.delay, G.mean_saccadeLatency, G.std_saccadeLatency, 'o-', 'LineWidth', 1.5);
% xlabel('Delay');
% ylabel('Mean Saccade Latency');
% prettyFig;

f = figure;
%histStyle_KKN(Tmdir.saccadeLatency)

%% cluster distributions

f = figure;
histStyle_KKN(Cleft.SPI)




%% checking that pursuit onset and saccade onset are aligned correctly

Tmdir = S.mdir1.tbl(S.mdir1.tbl.result=='CORRECT',:);
Tpurs = S.purs1.tbl(S.purs1.tbl.result=='CORRECT' & S.purs1.tbl.pursType=='pure',:);
Tpurs = Tpurs(isnan(Tpurs.msOffset) | Tpurs.msOffset<-100,:);

f1a = figure;
f1a.Position = [100 100 1500 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

% mdir, aligned to stim onset
nexttile
x = -300:1:500;
FR_WIN = [50 150];
fill([FR_WIN fliplr(FR_WIN)], [[1000 1000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')

for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyeVel{t}(1,:),Tmdir.eyeVel{t}(2,:));
    plot(x,rh(Tmdir.TARG_ON{t}(1)-300:Tmdir.TARG_ON{t}(1)+500),'k-')
end
xlabel('time aligned to stimulus onset (ms)')
ylabel('radial eye velocity (deg/s)')
ylim([0 1000])
title('visual epoch (MGS)')
prettyFig;

% purs, aligned to targ motion onset
nexttile
x = -300:1:500;
FR_WIN = [50 250];
%fill([FR_WIN fliplr(FR_WIN)], [[100 100] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tpurs)
    [~,rh] = cart2pol(Tpurs.eyeVel{t}(1,:),Tpurs.eyeVel{t}(2,:));
    [pursuit_onset, rxnTime, msOffset, csOnset, csVelocity, csPeak, csOffset, csAngle, csType] = detect_pursuitOnset(Tpurs.eyePos{t}, Tpurs.eyeVel{t}, Tpurs.PURSUIT_TARG_ON(t), S.purs1.params.crossingTime, Tpurs.pursuitSpeed(t), Tpurs.angle(t), 'CS_PREINT', 50, 'CS_POSTINT', 100, 'PLOT_TRACES', false);
    if isequal(csType,'pure')
        plot(x,rh(Tpurs.PURSUIT_TARG_ON(t)-300:Tpurs.PURSUIT_TARG_ON(t)+500),'k-')
    end
end
xlabel('time aligned to target motion onset (ms)')
ylabel('radial eye velocity (deg/s)')
ylim([0 50])
prettyFig;

% mdir, aligned to saccade onset
nexttile
x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[1000 1000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyeVel{t}(1,:),Tmdir.eyeVel{t}(2,:));
    plot(x,rh(Tmdir.SACCADE(t)-300:Tmdir.SACCADE(t)+500),'k-')
end
xlabel('time aligned to saccade onset (ms)')
ylabel('radial eye velocity (deg/s)')
ylim([0 1000])
title('motor epoch (MGS)')
prettyFig;

% purs, aligned to pursuit onset
nexttile
x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[100 100] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tpurs)
    [~,rh] = cart2pol(Tpurs.eyeVel{t}(1,:),Tpurs.eyeVel{t}(2,:));
    [pursuit_onset, rxnTime, msOffset, csOnset, csVelocity, csPeak, csOffset, csAngle, csType] = detect_pursuitOnset(Tpurs.eyePos{t}, Tpurs.eyeVel{t}, Tpurs.PURSUIT_TARG_ON(t), S.purs1.params.crossingTime, Tpurs.pursuitSpeed(t), Tpurs.angle(t), 'CS_PREINT', 50, 'CS_POSTINT', 100);
    if isequal(csType,'pure')
        plot(x,rh(Tpurs.pursuitOnset(t)-300:Tpurs.pursuitOnset(t)+500),'k-')
    end
    Tpurs.pursuitOnset(t) = pursuit_onset;
    Tpurs.pursuitLatency(t) = rxnTime;
    Tpurs.pursType(t) = categorical(string(csType));

end
xlabel('time aligned to pursuit onset (ms)')
ylabel('radial eye velocity (deg/s)')
ylim([0 50])
title('motor epoch (pursuit)')
prettyFig;

%%
f1a = figure;
f1a.Position = [100 100 1300 1000];
tl = tiledlayout(3,1);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

% % pos, stim onset
% nexttile
% x = -300:1:500;
% FR_WIN = [50 150];
% fill([FR_WIN fliplr(FR_WIN)], [[25 25] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
% hold on;
% xline(0,'k--')
% 
% for t=1:height(Tmdir)
%     [~,rh] = cart2pol(Tmdir.eyePos{t}(1,:),Tmdir.eyePos{t}(2,:));
%     plot(x,rh(Tmdir.TARG_ON{t}(1)-300:Tmdir.TARG_ON{t}(1)+500),'k-')
% end
% ylabel('radial eye position (deg)')
% xlim([-100 200])
% ylim([0 25])
% title('visual epoch (MGS)')
% prettyFig;

% pos, sacc onset
nexttile
x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[25 25] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyePos{t}(1,:),Tmdir.eyePos{t}(2,:));
    plot(x,rh(Tmdir.SACCADE(t)-300:Tmdir.SACCADE(t)+500),'k-')
end
xlim([-300 500])
ylim([0 25])
ylabel('radial eye position (deg)')
title('motor epoch (MGS)')
prettyFig;

% vel, stim onset
% nexttile
% x = -300:1:500;
% FR_WIN = [50 150];
% fill([FR_WIN fliplr(FR_WIN)], [[1000 1000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
% hold on;
% xline(0,'k--')
% 
% for t=1:height(Tmdir)
%     [~,rh] = cart2pol(Tmdir.eyeVel{t}(1,:),Tmdir.eyeVel{t}(2,:));
%     plot(x,rh(Tmdir.TARG_ON{t}(1)-300:Tmdir.TARG_ON{t}(1)+500),'k-')
% end
% ylabel('radial eye velocity (deg/sec)')
% xlim([-100 200])
% ylim([0 1000])
% prettyFig;

% vel, sacc onset
nexttile
x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[1000 1000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyeVel{t}(1,:),Tmdir.eyeVel{t}(2,:));
    plot(x,rh(Tmdir.SACCADE(t)-300:Tmdir.SACCADE(t)+500),'k-')
end
xlim([-300 500])
ylim([0 1000])
ylabel('radial eye velocity (deg/s)')
prettyFig;

% acc, stim onset
% nexttile
% x = -300:1:500;
% FR_WIN = [50 150];
% fill([FR_WIN fliplr(FR_WIN)], [[80000 80000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
% hold on;
% xline(0,'k--')
% 
% for t=1:height(Tmdir)
%     [~,rh] = cart2pol(Tmdir.eyeAcc{t}(1,:),Tmdir.eyeAcc{t}(2,:));
%     plot(x,rh(Tmdir.TARG_ON{t}(1)-300:Tmdir.TARG_ON{t}(1)+500),'k-')
% end
% xlabel('time aligned to stimulus onset (ms)')
% ylabel('radial eye acceleration (deg/sec2)')
% xlim([-100 200])
% ylim([0 80000])
% prettyFig;

% acc, sacc onset
nexttile
x = -300:1:500;
FR_WIN = [-50 50];
fill([FR_WIN fliplr(FR_WIN)], [[80000 80000] fliplr([0 0])], [130,130,130]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);
hold on;
xline(0,'k--')
for t=1:height(Tmdir)
    [~,rh] = cart2pol(Tmdir.eyeAcc{t}(1,:),Tmdir.eyeAcc{t}(2,:));
    plot(x,rh(Tmdir.SACCADE(t)-300:Tmdir.SACCADE(t)+500),'k-')
end
xlabel('time aligned to saccade onset (ms)')
ylabel('radial eye acceleration (deg/s^2)')
xlim([-300 500])
ylim([0 80000])
prettyFig;



%% plotting individual units rasters
S2 = S;
S2.purs1.tbl = Tpurs;

CLUSTER = 1;

%ia_mdirRasters(S2, 'ALIGN', 'stim', 'CLUSTER', CLUSTER)
%ia_mdirRasters(S2, 'ALIGN', 'sacc', 'CLUSTER', CLUSTER)
ia_mdirRasters(S2, 'ALIGN', 'fix_off', 'CLUSTER', CLUSTER)

%ia_pursRasters(S2, 'ALIGN', 'targ', 'PURE_ONLY', true, 'CLUSTER', CLUSTER)
ia_pursRasters(S2, 'ALIGN', 'purs', 'PURE_ONLY', true, 'CLUSTER', CLUSTER)




%% histograms

imec0 = S.kilosort(1).clusters;
imec1 = S.kilosort(2).clusters;

imec0 = imec0(imec0.mdir1_Hz>1 & imec0.purs1_Hz>1,:);
imec1 = imec1(imec1.mdir1_Hz>1 & imec1.purs1_Hz>1,:);

imec0 = imec0(imec0.vis_sel_dir>0.1 | imec0.sac_sel_dir>0.1 | imec0.pur_sel_dir>0.1,:);
imec1 = imec1(imec1.vis_sel_dir>0.1 | imec1.sac_sel_dir>0.1 | imec1.pur_sel_dir>0.1,:);

f1a = figure;
f1a.Position = [100 100 1200 400];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
values = imec0.purs1_Hz;
histStyle(values,'left','ypos','number of clusters',[min(values)-10 max(values)+10],[0 300], 1, 10, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.purs1_Hz;
histStyle(values,'left','ypos','number of clusters',[min(values)-10 max(values)+10],[0 300], 1, 10, [0 0 0]./255, [100 100 100]./255)

%%

f1a = figure;
f1a.Position = [100 100 1200 400];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
values = imec0.vis_sel_dir;
histStyle(values,'left','ypos','number of clusters',[min(values)-1 max(values)+1],[0 300], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.vis_sel_dir;
histStyle(values,'left','ypos','number of clusters',[min(values)-1 max(values)+1],[0 300], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

%%

f1a = figure;
f1a.Position = [100 100 1200 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
values = imec0.VMI;
histStyle(values,'left','VMI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.VMI;
histStyle(values,'right','VMI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec0.SPI;
histStyle(values,'left','SPI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.SPI;
histStyle(values,'right','SPI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)


%% 

f2a = figure;
f2a.Position = [100 100 1200 800];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
plot(imec0.VMI,imec0.SPI,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('SPI')
prettyFig;

nexttile
plot(imec1.VMI,imec1.SPI,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('SPI')
prettyFig;

%% 

f3a = figure;
f3a.Position = [100 100 1200 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
plot(imec0.SPI,imec0.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('SPI')
ylabel('ypos')
prettyFig;

nexttile
plot(imec1.SPI,imec1.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('SPI')
ylabel('ypos')
prettyFig;

nexttile
plot(imec0.VMI,imec0.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('ypos')
prettyFig;

nexttile
plot(imec1.VMI,imec1.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('ypos')
prettyFig;
