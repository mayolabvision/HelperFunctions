function ia_rfmpStim_CoverageMap(data, varargin)

p = inputParser;
addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
addParameter(p, 'TIME_BIN', [-400,0], @isnumeric);
addParameter(p, 'ALIGN_TO', 'TARG_ON', @ischar);
addParameter(p, 'ALIGN_IND', 1, @isnumeric);

parse(p, data, varargin{:});
data = p.Results.data;
TIME_BIN = p.Results.TIME_BIN;
ALIGN_TO = p.Results.ALIGN_TO;
ALIGN_IND = p.Results.ALIGN_IND;

%data = '/Volumes/lab_NHPdata-processed/kendra_scrappy_0174/tables/kendra_scrappy_0174.mat';

fprintf('\n------------------------------\n')
if ischar(data)
    [~, filename, ~] = fileparts(data);
    load(data,'S');
    fprintf(sprintf('\n----Data loaded for %s----\n',filename))
else
    S = data;
end


% Combine every task struct whose field name contains 'rfsa' (rfsa1, rfsa2, ...)
sFields = fieldnames(S);
rfsaFields = sFields(contains(sFields, 'rfsa', 'IgnoreCase', true));
if isempty(rfsaFields)
    error('No fields containing ''rfsa'' were found in the data struct.');
end
fprintf('Combining %d rfsa table(s): %s\n', numel(rfsaFields), strjoin(rfsaFields, ', '));

rfsaTbls = cell(numel(rfsaFields),1);
for i = 1:numel(rfsaFields)
    rfsaTbls{i} = S.(rfsaFields{i}).tbl;
end
thisTbl = vertcat(rfsaTbls{:});

thisTbl = thisTbl(thisTbl.result=='CORRECT',:);
nTrialsTotal = height(thisTbl);

% stimPos is stored in screen pixels -- convert to degrees of visual angle
screenDist = thisTbl.params(1,1).block.screenDistance;
pixPerCM = thisTbl.params(1,1).block.pixPerCM;
thisTbl.stimPos = cellfun(@(pos) pix2deg(pos, screenDist, pixPerCM), thisTbl.stimPos, 'UniformOutput', false);

% ALIGN_TO must be on the same ms clock as STIM_ON. It can either hold a
% single scalar per trial (e.g. SACCADE), or a cell array of multiple
% values per trial (e.g. TARG_ON, which can fire more than once) -- in the
% latter case ALIGN_IND picks which value within that trial's cell to align to
if ~ismember(ALIGN_TO, thisTbl.Properties.VariableNames)
    error('ALIGN_TO ''%s'' is not a column of thisTbl.', ALIGN_TO);
end
alignColAll = thisTbl.(ALIGN_TO);
if isnumeric(alignColAll)
    isCellAlign = false;
elseif iscell(alignColAll)
    isCellAlign = true;
else
    error('ALIGN_TO ''%s'' must be numeric-per-trial or a cell array of numeric vectors, not %s.', ALIGN_TO, class(alignColAll));
end

% Determine the full RF stimulus grid from the data itself, rather than
% assuming a fixed 12x12 layout
allPos = cell2mat(thisTbl.stimPos);
gridX = unique(allPos(:,1));
gridY = unique(allPos(:,2));
nX = numel(gridX);
nY = numel(gridY);
fprintf('Detected %d x %d stimulus grid (%d possible positions)\n', nX, nY, nX*nY);

angles = unique(thisTbl.angle);
nAngles = numel(angles);

% First pass: bin stimulus counts per angle so all subplots can share one
% color scale
countsByAngle = cell(nAngles,1);
nTrialsByAngle = zeros(nAngles,1);
for a = 1:nAngles
    angTbl = thisTbl(thisTbl.angle==angles(a),:);
    nTrialsByAngle(a) = height(angTbl);

    binPos = [];
    for t = 1:height(angTbl)
        pos = angTbl.stimPos{t};     % Nx2 [x y]
        onTimes = angTbl.STIM_ON{t}; % Nx1, ms relative to trial onset

        if isCellAlign
            alignVals = angTbl.(ALIGN_TO){t};
            if numel(alignVals) >= ALIGN_IND
                eventTime = alignVals(ALIGN_IND);
            else
                eventTime = NaN; % trial has fewer ALIGN_TO events than ALIGN_IND -- contributes nothing
            end
        else
            eventTime = angTbl.(ALIGN_TO)(t);
        end

        relOnTimes = onTimes - eventTime; % ms relative to alignEvent
        inBin = relOnTimes >= TIME_BIN(1) & relOnTimes < TIME_BIN(2);
        binPos = [binPos; pos(inBin,:)]; %#ok<AGROW>
    end

    [~, xi] = ismember(binPos(:,1), gridX);
    [~, yi] = ismember(binPos(:,2), gridY);
    valid = xi > 0 & yi > 0;
    countsByAngle{a} = accumarray([yi(valid), xi(valid)], 1, [nY, nX]);
end
cmax = max(cellfun(@(c) max(c(:)), countsByAngle));
cLim = [1, max(cmax, 1.0001)]; % 0-count squares are drawn white separately, so scale color to the nonzero counts

% Second pass: plot
figure('Name','RF Map Stimulus Coverage','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for a = 1:nAngles
    if isnumeric(angles)
        angleStr = sprintf('%g', angles(a));
    else
        angleStr = char(string(angles(a)));
    end

    counts = countsByAngle{a};
    countsPlot = counts;
    countsPlot(counts==0) = NaN; % draw 0-count squares as white, not the colormap's low end

    nexttile;
    h = imagesc(gridX, gridY, countsPlot, cLim);
    set(h, 'AlphaData', ~isnan(countsPlot));
    set(gca,'YDir','normal','Color','w');
    axis square;
    colorbar;
    colormap(gca, 'parula');
    xlabel('X Position (deg)');
    ylabel('Y Position (deg)');
    title(sprintf('Angle = %s (n = %d trials)', angleStr, nTrialsByAngle(a)));
    prettyFig;
end

if isCellAlign
    alignStr = sprintf('%s(%d)', ALIGN_TO, ALIGN_IND);
else
    alignStr = ALIGN_TO;
end
sgtitle(sprintf('RF Stimulus Coverage | Time Bin [%d, %d] ms aligned to %s | Total Correct Trials = %d', ...
    TIME_BIN(1), TIME_BIN(2), alignStr, nTrialsTotal), 'Interpreter', 'none');

