function raster_sdf(sptimes,timewindow,sigma,line_color,sem_shade)

% raster for a single neuron (each row is a trial, each col is a spike time) 

%%%%%%%%%% Inputs: %%%%%%%%%%%%%
% sptimes --> 1 x N cell array (N = number of trials)
   % sptimes{n} = spike times in specified window
% timewindow --> time points (not indices) to include (e.g. [-100 300])
% sigma --> width of gaussian/window [ms] (e.g. 5)
% line_color --> color of lines (e.g. [0 0 0]./255)
% shade_color --> shade of lines (e.g. [153 153 153]./255)

if nargin < 4
    line_color = [0 0 0]./255;
    sem_shade = [153 153 153]./255;
end

% f = figure; %figure('Units','normalized','Position',[0 0 .3 1])
% f.Position = [100 100 900 600];
ax = tiledlayout(2,1);
ax.TileSpacing = 'compact';
ax.Padding = 'compact';

nexttile(ax)
% For all trials...
for iTrial = 1:length(sptimes)
                  
    spks            = sptimes{iTrial}';         % Get all spikes of respective trial 
%     if size(spks,1)==1
%         spks = spks';
%     end
    xspikes         = repmat(spks,3,1);         % Replicate array
    yspikes      	= nan(size(xspikes));       % NaN array
    
    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                % Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end
    
    plot(xspikes, yspikes, 'Color', line_color , 'LineWidth', 2)
    hold on
end
xline(0,'k--','linewidth',1)

xlim(timewindow)
ylim([0 length(sptimes)])
yticks([0 length(sptimes)])
%xticks([min(timewindow) 0 max(timewindow)])

%ylabel('trials')
%title(ax,[sprintf('%d',condition) char(176)],'fontweight','bold','fontsize',14)

prettyFig;

%% Spike density function

tstep     	= 1;                                          % Resolution for SDF [ms]                   
time     	= tstep+timewindow(1):tstep:timewindow(2);    % Time vector

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
x = ((1:size(sdf,2)) + timewindow(1));
fill([x fliplr(x)], [yu fliplr(yl)], sem_shade, 'linestyle','none','FaceAlpha',0.5)
hold on
plot(x,mn, 'Color', line_color, 'LineWidth', 1.5)
xline(0,'k--','linewidth',1)
mVal = max(yu) + round(max(yu)*.1);

xlim(timewindow)
%ylim([0 mVal])
ylim([0 60])
%xticks([min(timewindow) max(timewindow)])

%ylabel('FR [Hz]')
prettyFig;

end