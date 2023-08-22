clear
clc

addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions'))

datafolder = '/Users/kendranoneman/Projects/mayo/local_data/rig_setup';
figfolder = '/Users/kendranoneman/Figures/patrick_r01';
datafile = 'sb05pursA65650033';
kernel = 50;

tbl = extract_eyeTraces(datafolder,datafile,kernel);

%% Detect pursuit onset, pursuit offset, catch-up saccades

data = tbl(tbl.trialOutcome=="CORRECT",:);
dirs = sort(unique(data.angle));

preint = 350;
x1 = (1:preint+1);
x2 = (1:preint+1) - preint;

f1a = figure;
f1a.Position = [100 100 1400 900];
tl = tiledlayout(3,2);
tl.TileSpacing = 'loose'; 
tl.Padding = 'compact';

one_dir = data; %(data.angle==dirs(d),:);
[~,Rbeg] = cellfun(@(q,r) cart2pol(q(1,r+100:r+350),q(2,r+100:r+350)), one_dir.eyeAcc,num2cell(one_dir.TARG_ON), 'uni', 0);
[~,Rend] = cellfun(@(q,r) cart2pol(q(1,r-350:r),q(2,r-350:r)), one_dir.eyeAcc,num2cell(one_dir.REWARD), 'uni', 0);

saccBeg = cellfun(@(q) sum(q>=1000,2)>0, Rbeg, 'uni', 1);
saccEnd = cellfun(@(q) sum(q>=1000,2)>0, Rend, 'uni', 1);

% POSITION
a(1) = nexttile;
for t=1:height(one_dir)
    [~,R] = cart2pol(one_dir.eyePos{t}(1,one_dir.TARG_ON(t):one_dir.TARG_ON(t)+preint),one_dir.eyePos{t}(2,one_dir.TARG_ON(t):one_dir.TARG_ON(t)+preint));

    if saccBeg(t)==0
        plot(x1,R,'k-')
    %else
        %plot(x1,R,'b-')
    end
    hold on
end
xlim([100 preint])
xline(100,'k--')
xline(200,'k--')
title([sprintf('%2.0f',(sum(saccBeg)/length(saccBeg))*100),'% trials'])
prettyFig;

aa(1) = nexttile;
for t=1:height(one_dir)
    [~,R] = cart2pol(one_dir.eyePos{t}(1,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),one_dir.eyePos{t}(2,one_dir.REWARD(t)-preint:one_dir.REWARD(t)));

    if saccEnd(t)==0
        plot(x2,R,'k-')
%     else
%         plot(x2,R,'b-')
    end
    hold on
end
xlim([-preint -100])
xline(-300,'k--')
xline(100,'k--')
title([sprintf('%2.0f',(sum(saccEnd)/length(saccEnd))*100),'% trials'])
prettyFig;

% VELOCITY
a(2) = nexttile;
for t=1:height(one_dir)
    [~,R] = cart2pol(smoothdata(one_dir.eyeVel{t}(1,one_dir.TARG_ON(t):one_dir.TARG_ON(t)+preint),'gaussian',50),smoothdata(one_dir.eyeVel{t}(2,one_dir.TARG_ON(t):one_dir.TARG_ON(t)+preint),'gaussian',50));

    if saccBeg(t)==0
        plot(x1,R,'k-')
%     else
%         plot(x1,R,'b-')
    end
    hold on
end
xlim([100 preint])
ylim([0 20])
xline(100,'k--')
xline(200,'k--')
prettyFig;

aa(2) = nexttile;
for t=1:height(one_dir)
    [~,R] = cart2pol(smoothdata(one_dir.eyeVel{t}(1,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),'gaussian',50),smoothdata(one_dir.eyeVel{t}(2,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),'gaussian',50));

    if saccEnd(t)==0
        plot(x2,R,'k-')
%     else
%         plot(x2,R,'b-')
    end
    hold on
end
xlim([-preint -100])
ylim([0 20])
xline(-300,'k--')
xline(100,'k--')
prettyFig;

% ACCELERATION
a(3) = nexttile;
for t=1:height(one_dir)
    [~,R] = cart2pol(one_dir.eyeAcc{t}(1,one_dir.TARG_ON(t):one_dir.TARG_ON(t)+preint),one_dir.eyeAcc{t}(2,one_dir.TARG_ON(t):one_dir.TARG_ON(t)+preint));

    if saccBeg(t)==0
        plot(x1,R,'k-')
%     else
%         plot(x1,R,'b-')
    end
    hold on
end
xlim([100 preint])
xline(100,'k--')
xline(200,'k--')
xlabel('time aligned to target motion onset (ms)')
prettyFig;

aa(3) = nexttile;
for t=1:height(one_dir)
    [~,R] = cart2pol(one_dir.eyeAcc{t}(1,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),one_dir.eyeAcc{t}(2,one_dir.REWARD(t)-preint:one_dir.REWARD(t)));

    if saccEnd(t)==0
        plot(x2,R,'k-')
%     else
%         plot(x2,R,'b-')
    end
    hold on
end
xlim([-preint -100])
xline(-300,'k--')
xline(100,'k--')
xlabel('time aligned to target offset (ms)')
prettyFig;

linkaxes(a,'x')
linkaxes(aa,'x')
title(tl,sprintf('8 directions (%d trials)',length(saccBeg)),'fontsize',16)

savebigPDF(f1a, sprintf('%s/alldirs_catchupSaccades_justPurePursuit_smoothed.pdf',figfolder))


%%
data = tbl(tbl.trialOutcome=="CORRECT",:);
dirs = sort(unique(data.angle));

one_dir = data; %(data.angle==dirs(d),:);

preints = [450,425,400,375,350,317,300];

f1a = figure;
f1a.Position = [100 100 1400 1400];
tl = tiledlayout(length(preints),2);
tl.TileSpacing = 'loose'; 
tl.Padding = 'compact';

for p=1:length(preints)
    preint = preints(p);
    x2 = (1:preint+1) - preint;
    postint = 0;
    
    [~,Rend] = cellfun(@(q,r) cart2pol(q(1,r-preint:r+postint),q(2,r-preint:r+postint)), one_dir.eyeAcc,num2cell(one_dir.REWARD), 'uni', 0);
    saccEnd = cellfun(@(q) sum(q>=1000,2)>0, Rend, 'uni', 1);
    
    nexttile
    for t=1:height(one_dir)
        [~,R] = cart2pol(smoothdata(one_dir.eyeVel{t}(1,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),'gaussian',50),smoothdata(one_dir.eyeVel{t}(2,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),'gaussian',50));
    
        if saccEnd(t)==0
            plot(x2,R,'k-')
        else
            plot(x2,R,'b-')
        end
        hold on
    end
    xlim([-preint 0])
    title([sprintf('%2.0f',(sum(saccEnd)/length(saccEnd))*100), '%'],'fontsize',10)
    prettyFig;
    
    nexttile
    for t=1:height(one_dir)
        [~,R] = cart2pol(smoothdata(one_dir.eyeVel{t}(1,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),'gaussian',50),smoothdata(one_dir.eyeVel{t}(2,one_dir.REWARD(t)-preint:one_dir.REWARD(t)),'gaussian',50));
    
        if saccEnd(t)==0
            plot(x2,R,'k-')
    %     else
    %         plot(x2,R,'b-')
        end
        hold on
    end
    xlim([-preint 0])
    title([sprintf('%2.0f',length(saccEnd)-sum(saccEnd)), ' pure trials'],'fontsize',10)
    ylim([0 20])
    prettyFig;
end
%title(tl,sprintf('8 directions (%d trials)',length(saccBeg)),'fontsize',16)

savebigPDF(f1a, sprintf('%s/alldirs_catchupSaccades_differentWindowSizes.pdf',figfolder))

%% Simplifying data struct
preint = 200; % how many seconds before targ onset
postint = 1400; % how many ms before reward 

data = dat(cellfun(@isnan, {dat.result}.', 'uni', 1));

conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), {data.text}.', 'uni', 0);
conditions = vertcat(conditions{:}); conditions = conditions(:,2:5);

eye = {data.eyedata}.'; 
times = cellfun(@(q) round(q(1)*1000:q(2)*1000), {data.time}.', 'uni', 0);
trialcodes = {data.trialcodes}.';

inds = cellfun(@(q,r) find(q == round(r(r(:,2)==70,3)*1000)), times, trialcodes, 'uni', 0);
eyeTraces = cellfun(@(x,y) smoothdata(x.trial(1:2,y-preint:y+postint),2,'gaussian',60), eye, inds, 'uni', 0);

HEPos = cellfun(@(q) q(1,:), eyeTraces, 'uni', 0);

x = 1:length(HEPos{1});
HEVel = cellfun(@(q) smoothdata((gradient(q(:)) ./ gradient(x(:)./1000))','gaussian',20), HEPos, 'uni', 0);
HEAcc = cellfun(@(q) smoothdata((gradient(q(:)) ./ gradient(x(:)./1000))','gaussian',20), HEVel, 'uni', 0);

tbl = cell2table([conditions cellfun(@(x) {x}, HEPos, 'uni', 0) cellfun(@(x) {x}, HEVel, 'uni', 0) cellfun(@(x) {x}, HEAcc, 'uni', 0)],'VariableNames',["fixDuration","pursuitSpeed","angle","jumpSize","HE_eyePosition","HE_eyeVelocity","HE_eyeAcceleration"]);

%% FULL TRACES
angles = flip(unique(tbl.angle));

f1a = figure;
f1a.Position = [100 100 1400 900];
tl = tiledlayout(3,length(angles));
tl.TileSpacing = 'loose'; 
tl.Padding = 'compact';

x = (1:length(HEPos{1})) - preint;

for b=1:3
    if b==1
        dd = tbl.HE_eyePosition;
        ylab = 'position (deg)';
    elseif b==2
        dd = tbl.HE_eyeVelocity;
        ylab = 'velocity (deg/s)';
    else
        dd = tbl.HE_eyeAcceleration;
        ylab = 'acceleration (deg/s^2)';
    end
    ax = [];
    for a=1:length(angles)
        ax(a) = nexttile;
        trc = dd(tbl.angle==angles(a),:); trc = vertcat(trc{:});
    
        plot(x,trc,'k-','linewidth',1)
        xlim([-preint,postint])
        if b==1
            title(sprintf('%d deg (%d total trials)',angles(a),size(dd(tbl.angle==angles(a),:),1)))
        end
        if a==1
            ylabel(ylab)
        end
        prettyFig;
        linkaxes(ax,'xy')
    end
end

xlabel(tl,'time aligned to target onset (ms)','fontsize',22)
title(tl,datafile,'fontsize',24)
saveas(f1a,sprintf('%s/%s_full.png',datafolder,datafile))

%% ZOOMED IN TRACES
f1b = figure;
f1b.Position = [100 100 1400 900];
tl = tiledlayout(3,2);
tl.TileSpacing = 'loose'; 
tl.Padding = 'compact';

x = (1:length(HEPos{1})) - preint;
window = [100 200]+preint;

acc = tbl.HE_eyeAcceleration; acc = vertcat(acc{:});
sacc = sum(abs(acc(:,window(1):window(2))) > 1000,2) > 0;

window = [50 300];
for b=1:3
    if b==1
        dd = tbl.HE_eyePosition;
        ylab = 'position (deg)';
    elseif b==2
        dd = tbl.HE_eyeVelocity;
        ylab = 'velocity (deg/s)';
    else
        dd = tbl.HE_eyeAcceleration;
        ylab = 'acceleration (deg/s^2)';
    end
    dd = vertcat(dd{:});

    ax = [];
    for a=1:length(angles)
        ax(a) = nexttile;
        plot(x,dd(tbl.angle==angles(a) & ~sacc,:),'k-','linewidth',1);
        hold on
        plot(x,dd(tbl.angle==angles(a) & sacc,:),'b-','linewidth',1);
        xlim(window)
        if b==1
            percSacc = (size(dd(tbl.angle==angles(a) & sacc,:),1)/size(dd(tbl.angle==angles(a),:),1))*100;
            title(sprintf('%d deg (%d total trials)',angles(a),size(dd(tbl.angle==angles(a),:),1)))
            subtitle([sprintf('(%0.2f',percSacc),'% trials w/ saccades)'])
        end
        if a==1
            ylabel(ylab)
        end
        prettyFig;
        linkaxes(ax,'xy')
    end
end

xlabel(tl,'time aligned to target onset (ms)','fontsize',22)
title(tl,datafile,'fontsize',24)
saveas(f1b,sprintf('%s/%s_zoom.png',datafolder,datafile))
