function raster_sdf(spike_times, varargin)
% RASTER_SDF plots a raster and spike density function (SDF) for a single neuron.
%
% Inputs:
%   sptimes      - 1xN cell array (N = number of trials), each cell contains spike times
%   TIME_WINDOW   - Time points (not indices) to include (e.g. [-100 300])
%   SIGMA        - Width of Gaussian/window [ms] (e.g. 5)
%
% Optional Inputs (via varargin):
%   line_color   - Color of the SDF line (default: [0 0 255]./255)
%   sem_shade    - Color of the SEM shading (default: [178 178 255]./255)
%   rast_shade   - Color of raster plot lines (default: [166 166 166]./255)


% Parse varargin inputs
p = inputParser;
addRequired(p, 'spike_times', @iscell);
addOptional(p,'TIME_WINDOW', [-300,500], @(x) isnumeric(x) && length(x) == 2) % ms
addOptional(p,'SIGMA', 10, @isnumeric)
addOptional(p,'LINE_COLOR', [0 0 255]./255, @(x) (isnumeric(x) && length(x) == 3) || (iscell(x) && length(x) == length(spike_times)));
addOptional(p,'SEM_SHADE', [178 178 255]./255, @(x) (isnumeric(x) && length(x) == 3) || (iscell(x) && length(x) == length(spike_times)));
addOptional(p,'TICK_COLOR', [30 30 30]./255, @(x) (isnumeric(x) && length(x) == 3) || (iscell(x) && length(x) == length(spike_times)));
addOptional(p,'TICK_LENGTH', ceil(numel(spike_times)/50), @isnumeric)
addOptional(p,'TICK_WIDTH', 2, @isnumeric)
addOptional(p,'FR_WINDOW', [], @isnumeric)

p.parse(spike_times, varargin{:});
spike_times = p.Results.spike_times;
TIME_WINDOW = p.Results.TIME_WINDOW;
SIGMA = p.Results.SIGMA;
LINE_COLOR = p.Results.LINE_COLOR;
SEM_SHADE  = p.Results.SEM_SHADE;
TICK_COLOR = p.Results.TICK_COLOR;
TICK_LENGTH = p.Results.TICK_LENGTH;
TICK_WIDTH = p.Results.TICK_WIDTH;
FR_WINDOW = p.Results.FR_WINDOW;

if iscell(LINE_COLOR)
    unique_line_colors = num2cell(unique(cell2mat(LINE_COLOR),'rows'),2);
    unique_sem_shades = num2cell(unique(cell2mat(SEM_SHADE),'rows'),2);
end

%% Raster Plot
yyaxis right; 
ax = gca;
ax.YColor = 'k';

yl_back = [0 0]; yu_back = [numel(spike_times) numel(spike_times)];
fill([FR_WINDOW fliplr(FR_WINDOW)], [yu_back fliplr(yl_back)], [230,230,230]./255, 'linestyle', 'none', 'FaceAlpha', 0.25);

hold on;

for iTrial = 1:length(spike_times)
    spks = spike_times{iTrial};
    if size(spks,2)==1
        spks = spks';
    end
    xspikes = repmat(spks, 3, 1);
    yspikes = nan(size(xspikes));
    
    %if ~isempty(yspikes)
    yspikes(1, :) = iTrial - TICK_LENGTH;
    yspikes(2, :) = iTrial;
    
    %end
    
    if ~iscell(TICK_COLOR)
        plot(xspikes, yspikes, '-', 'Color', TICK_COLOR, 'LineWidth', TICK_WIDTH);
    else
        plot(xspikes, yspikes, '-', 'Color', TICK_COLOR{iTrial}, 'LineWidth', TICK_WIDTH);
    end
    hold on;
end

xline(0, 'k-', 'linewidth', 1);
xlim(TIME_WINDOW);
ylim([0 length(spike_times)]);
if ~iscell(TICK_COLOR)
    yticks([0 length(spike_times)]);
else
    block_changes = find([false; ~cellfun(@isequal, TICK_COLOR(1:end-1), TICK_COLOR(2:end))])';
    if numel(block_changes) > 1
        yticks([0 length(spike_times)]);
    else
        yticks([0 block_changes-1 length(spike_times)]);
        plot([TIME_WINDOW(2)+5 TIME_WINDOW(2)+5], [1 block_changes-1], '-', 'Color', unique_line_colors{1}, 'LineWidth', 2)
        plot([TIME_WINDOW(2)+5 TIME_WINDOW(2)+5], [block_changes length(spike_times)], '-', 'Color', unique_line_colors{2}, 'LineWidth', 2)
        xlim([TIME_WINDOW(1) TIME_WINDOW(2)+10])
    end
end
ylabel('trials','Rotation',270);
prettyFig;

%% Spike Density Function
yyaxis left;
ax = gca;
ax.YColor = 'k';

% Time vector
tstep = 1;
time = tstep + TIME_WINDOW(1) : tstep : TIME_WINDOW(2);

sdf = zeros(length(spike_times), length(time));
for iTrial = 1:length(spike_times)
    spks = spike_times{iTrial};
    if size(spks,2)==1
        spks = spks';
    end
    if isempty(spks)
        continue;
    end
    
    gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / SIGMA) .^ 2) ./ (SIGMA * sqrt(2 * pi)), spks, 'UniformOutput', false);
    sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
end

% Compute mean and SEM
x = (1:size(sdf, 2)) + TIME_WINDOW(1);

if ~iscell(LINE_COLOR)
    [mn, ~, yu, yl] = sem_errorbar(sdf .* 1000);
    fill([x fliplr(x)], [yu fliplr(yl)], SEM_SHADE, 'linestyle', 'none', 'FaceAlpha', 0.5);
    hold on;
    plot(x, mn, '-', 'Color', LINE_COLOR, 'LineWidth', 2);
else
    for dd = 1:numel(unique_line_colors)
        sdf2 = sdf(cellfun(@(q) isequal(q, unique_line_colors{dd}), LINE_COLOR),:);
        [mn, ~, yu, yl] = sem_errorbar(sdf2 .* 1000);
        fill([x fliplr(x)], [yu fliplr(yl)], unique_sem_shades{dd}, 'linestyle', 'none', 'FaceAlpha', 0.5);
        hold on;
        plot(x, mn, '-', 'Color', unique_line_colors{dd}, 'LineWidth', 2);


    end
end

ylabel('firing rate (Hz)');
prettyFig;

end
