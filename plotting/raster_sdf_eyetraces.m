function [ax1,ax2] = raster_sdf_eyetraces(sptimes,eyevels,timewindow,align,sigma,condition,p,col,sha)
% 
% Inputs:
% sptimes --> 1 x D cell array (D = number of neurons)
   % sptimes{d} = N x 1 double array (N = number of trials)
% timewindow --> time before/after "onset" in ms (e.g. [-100 300])
% align --> "stim" (target motion onset), "pursuit" (pursuit onset), or "saccade" (catch-up saccade onset)

% f = figure; %figure('Units','normalized','Position',[0 0 .3 1])
% f.Position = [100 100 900 600];
ax = tiledlayout(p,3,1);
ax.TileSpacing = 'compact';
ax.Padding = 'compact';

nexttile(ax)
% For all trials...
for iTrial = 1:length(sptimes)
                  
    spks            = sptimes{iTrial}';         % Get all spikes of respective trial    
    xspikes         = repmat(spks,3,1);         % Replicate array
    yspikes      	= nan(size(xspikes));       % NaN array
    
    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                % Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end
    
    plot(xspikes, yspikes, 'Color', col, 'LineWidth', 2)
    hold on
end
xline(0,'k--','linewidth',1)

xlim(timewindow)
ylim([0 length(sptimes)])
yticks([0 length(sptimes)])
%xticks([min(timewindow) 0 max(timewindow)])

ylabel('trials')
title(ax,[sprintf('%d',condition) char(176)],'fontweight','bold','fontsize',14)

prettyFig;

%% Spike density function

tstep     	= 1;                     % Resolution for SDF [ms]
%sigma   	= sigma;                     % Width of gaussian/window [ms]
time     	= tstep+timewindow(1):tstep:timewindow(2);        % Time vector

for iTrial = 1:length(sptimes)
    spks    = []; 
    gauss   = []; 
    spks   	= sptimes{iTrial}';          % Get all spikes of respective trial
    
    if isempty(spks)            
        out	= zeros(1,length(time));    % Add zero vector if no spikes
    else
        
        % For every spike
        for iSpk = 1:length(spks)
            
            % Center gaussian at spike time
            mu              = spks(iSpk);
            
            % Calculate gaussian
            p1              = -.5 * ((time - mu)/sigma) .^ 2;
            p2              = (sigma * sqrt(2*pi));
            gauss(iSpk,:)   = exp(p1) ./ p2;
            
        end
        
        % Sum over all distributions to get spike density function
        sdf(iTrial,:)       = sum(gauss,1);
    end
end

[mn,~,yu,yl] = sem_errorbar(sdf.*1000);

% Average response
ax1 = nexttile(ax);
x = (1:size(sdf,2)) + timewindow(1);
fill([x fliplr(x)], [yu fliplr(yl)], sha, 'linestyle','none','FaceAlpha',0.5)
hold on
plot(x,mn, 'Color', col, 'LineWidth', 1.5)
xline(0,'k--','linewidth',1)
mVal = max(yu) + round(max(yu)*.1);

xlim(timewindow)
ylim([0 mVal])
%xticks([min(timewindow) max(timewindow)])

% if isequal(align,'stim')
%     xlabel('Time aligned to target motion onset [ms]')
% elseif isequal(align,'pursuit')
%     xlabel('Time aligned to pursuit onset [ms]')
% elseif isequal(align,'saccade')
%     xlabel('Time aligned to saccade onset [ms]')
% else
%     xlabel('Time [ms]')
% end
ylabel('FR [Hz]')

prettyFig;

%% Eye traces

[mn,~,yu,yl] = sem_errorbar(eyevels);

% Average response
ax2 = nexttile(ax);
x = (1:length(eyevels)) - 300;
fill([x fliplr(x)], [yu fliplr(yl)], sha, 'linestyle','none','FaceAlpha',0.5)
hold on
plot(x,mn, 'Color', col, 'LineWidth', 1.5)
xline(0,'k--','linewidth',1)
%mVal = max(mn) + round(max(mn)*.1);

xlim(timewindow)
%ylim([0 mVal])
%xticks([min(timewindow) max(timewindow)])
yticks([0 15])

if isequal(align,'stim')
    xlabel('Time aligned to target motion onset [ms]')
elseif isequal(align,'pursuit')
    xlabel('Time aligned to pursuit onset [ms]')
elseif isequal(align,'saccade')
    xlabel('Time aligned to saccade onset [ms]')
else
    xlabel('Time [ms]')
end
ylabel('eye velocity [deg/s]')

prettyFig;


end