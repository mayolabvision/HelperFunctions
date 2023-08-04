clear
clc

addpath(genpath('/Users/kendranoneman/Projects/mayo/packages/nevutils'))
datafolder = '/Users/kendranoneman/Projects/mayo/HelperFxns/process_recordings/example_data';

datafile = 'sb01pursA65650026';
[dat,hdr] = nev2dat(sprintf('%s/%s',datafolder,datafile),'readNS5',true,'convertEyes',true);

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
