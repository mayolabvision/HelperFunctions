function raster_sdf(sptimes, timewindow, sigma, varargin)
% RASTER_SDF plots a raster and spike density function (SDF) for a single neuron.
%
% Inputs:
%   sptimes      - 1xN cell array (N = number of trials), each cell contains spike times
%   timewindow   - Time points (not indices) to include (e.g. [-100 300])
%   sigma        - Width of Gaussian/window [ms] (e.g. 5)
%
% Optional Inputs (via varargin):
%   line_color   - Color of the SDF line (default: [0 0 255]./255)
%   sem_shade    - Color of the SEM shading (default: [178 178 255]./255)
%   rast_shade   - Color of raster plot lines (default: [166 166 166]./255)

% Set default values for optional parameters
default_line_color = [0 0 255] ./ 255;
default_sem_shade  = [178 178 255] ./ 255;
default_rast_shade = [10 10 10] ./ 255;

% Parse varargin inputs
p = inputParser;
p.addOptional('line_color', default_line_color, @(x) isnumeric(x) && length(x) == 3);
p.addOptional('sem_shade', default_sem_shade, @(x) isnumeric(x) && length(x) == 3);
p.addOptional('rast_shade', default_rast_shade, @(x) isnumeric(x) && length(x) == 3);
p.parse(varargin{:});
line_color = p.Results.line_color;
sem_shade  = p.Results.sem_shade;
rast_shade = p.Results.rast_shade;

%% Raster Plot
yyaxis right; 
ax = gca;
ax.YColor = 'k';

for iTrial = 1:length(sptimes)
    spks = sptimes{iTrial};
    if size(spks,2)==1
        spks = spks';
    end
    xspikes = repmat(spks, 3, 1);
    yspikes = nan(size(xspikes));
    
    if ~isempty(yspikes)
        yspikes(1, :) = iTrial - 1;
        yspikes(2, :) = iTrial;
    end
    
    plot(xspikes, yspikes, '-', 'Color', rast_shade, 'LineWidth', 1);
    hold on;
end

xline(0, 'k-', 'linewidth', 2);
xlim(timewindow);
ylim([0 length(sptimes)]);
yticks([0 length(sptimes)]);
ylabel('trials','Rotation',270);
prettyFig;

%% Spike Density Function
yyaxis left;
ax = gca;
ax.YColor = 'k';

% Time vector
tstep = 1;
time = tstep + timewindow(1) : tstep : timewindow(2);

sdf = zeros(length(sptimes), length(time));
for iTrial = 1:length(sptimes)
    spks = sptimes{iTrial};
    if size(spks,2)==1
        spks = spks';
    end
    if isempty(spks)
        continue;
    end
    
    gauss = arrayfun(@(mu) exp(-.5 * ((time - mu) / sigma) .^ 2) ./ (sigma * sqrt(2 * pi)), spks, 'UniformOutput', false);
    sdf(iTrial, :) = sum(cell2mat(gauss'), 1);
end

% Compute mean and SEM
x = (1:size(sdf, 2)) + timewindow(1);
[mn, ~, yu, yl] = sem_errorbar(sdf .* 1000);
fill([x fliplr(x)], [yu fliplr(yl)], sem_shade, 'linestyle', 'none', 'FaceAlpha', 0.5);
hold on;
plot(x, mn, '-', 'Color', line_color, 'LineWidth', 2);

ylabel('firing rate (Hz)');
prettyFig;

end
