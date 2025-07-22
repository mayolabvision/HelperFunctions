%%%%%%%%%% Code to plot polar plot of MGS eye positions %%%%%%%%%%
% For Ani's GRC presentation
% S struct contains data table w/ task codes and behavior

%% 1. Load in processed table for example session 
table_path = '/Volumes/SHARED_STUFF/tables_behavior/kendra_scrappy_0124a_g0_behavOnly';
load(table_path,'S');

%% 2. Plot eye position for each correct trial in polar coordinates
% Duration of time (in ms) to include PRE and POST initiation of saccade
% (using SACCADE code in Ex, which is when his eyes left the fix window)
PREINT = 10; POSTINT = 75;

% Use only correct trials for plotting
T = S.mdir1.tbl(S.mdir1.tbl.result=='CORRECT',:);

f1 = figure;

% Loop through each trial in the table
for t = 1:height(T)
    this_angle    =  T.angle(t);   % angle in degrees
    this_saccTime =  T.SACCADE(t); % time (ms) of saccade

    % Pull out eye traces for trial and convert to polar coordinates
    this_eyePos = T.eyePos{t}; % HE, VE eye traces
    [theta,rho] = cart2pol(T.eyePos{t}(1,:),T.eyePos{t}(2,:));

    % Adding a check here to not include trials where his first saccade was in the wrong direction or way beyond the bounds of the screen
    if max(rho(this_saccTime-PREINT:this_saccTime+POSTINT)) < 40 &&  abs(wrapTo180(mean(rad2deg(theta(this_saccTime-PREINT:this_saccTime+POSTINT)))-this_angle)) < 45
        polarplot(theta(this_saccTime-PREINT:this_saccTime+POSTINT),rho(this_saccTime-PREINT:this_saccTime+POSTINT), ...
            'color',radial_colormap(this_angle, 'SHIFT_DEG', 90),'linewidth',1)
        hold on;
    end
    prettyFig;
end

