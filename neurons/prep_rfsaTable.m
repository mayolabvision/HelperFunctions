function [thisTbl, alignStr, nTrialsTotal, nStimuli] = prep_rfsaTable(S, ALIGN_TO, ALIGN_IND, TIME_BIN, ANGLE)
% prep_rfsaTable — Shared data prep for the ia_rfmpStim_RFMap* family.
% Combines every task struct whose field name contains 'rfsa', filters to
% CORRECT trials (and to ANGLE if specified), converts stimPos from screen
% pixels to degrees of visual angle, and restricts each trial's stimPos/
% STIM_ON to only the RF flashes that fall within TIME_BIN relative to
% ALIGN_TO(ALIGN_IND) -- flashes elsewhere in the trial (e.g. after the
% target appears) are dropped before anything else uses this table.
%
% INPUTS:
%   S         : session struct (already loaded, not a file path)
%   ALIGN_TO  : column in thisTbl to align STIM_ON to -- either a single
%               numeric time per trial (e.g. SACCADE), or a cell array of
%               multiple values per trial (e.g. TARG_ON)
%   ALIGN_IND : if ALIGN_TO is a cell column, which value within that
%               trial's cell to align to
%   TIME_BIN  : [start end] ms relative to ALIGN_TO(ALIGN_IND) -- only
%               STIM_ON flashes in this window are kept
%   ANGLE     : if not NaN, only use trials where thisTbl.angle == ANGLE
%
% OUTPUTS:
%   thisTbl      : the prepped table, with 'conditions' built (the nested-
%                  cell format format_tableToRFMap expects) and stimPos/
%                  STIM_ON restricted to the qualifying flashes
%   alignStr     : display string for the alignment event, e.g. 'TARG_ON(1)'
%   nTrialsTotal : number of correct (+ANGLE-filtered) trials before the
%                  TIME_BIN restriction
%   nStimuli     : total qualifying RF flashes across all trials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if ~isnan(ANGLE)
    if ~ismember(ANGLE, thisTbl.angle)
        error('ANGLE %g not found among thisTbl.angle values: %s', ANGLE, mat2str(unique(thisTbl.angle)'));
    end
    thisTbl = thisTbl(thisTbl.angle==ANGLE,:);
    fprintf('Restricting to angle == %g (%d trials)\n', ANGLE, height(thisTbl));
end
nTrialsTotal = height(thisTbl);

% stimPos is stored in screen pixels -- convert to degrees of visual angle
screenDist = thisTbl.params(1,1).block.screenDistance;
pixPerCM = thisTbl.params(1,1).block.pixPerCM;
thisTbl.stimPos = cellfun(@(pos) pix2deg(pos, screenDist, pixPerCM), thisTbl.stimPos, 'UniformOutput', false);

% ALIGN_TO must be on the same ms clock as STIM_ON
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

% Keep only the RF flashes (stimPos/STIM_ON) that fall within TIME_BIN
% relative to ALIGN_TO(ALIGN_IND). format_tableToRFMap expects a
% 'conditions' column: per trial, a 1xN cell array where each cell is a
% 1x2 [x y] stimulus position (paired by index with that trial's STIM_ON).
thisTbl.conditions = cell(height(thisTbl),1);
for t = 1:height(thisTbl)
    pos = thisTbl.stimPos{t};     % Nx2 [x y]
    onTimes = thisTbl.STIM_ON{t}; % Nx1, ms relative to trial onset

    if isCellAlign
        alignVals = thisTbl.(ALIGN_TO){t};
        if numel(alignVals) >= ALIGN_IND
            eventTime = alignVals(ALIGN_IND);
        else
            eventTime = NaN; % trial has fewer ALIGN_TO events than ALIGN_IND -- contributes nothing
        end
    else
        eventTime = thisTbl.(ALIGN_TO)(t);
    end

    relOnTimes = onTimes - eventTime;
    inBin = relOnTimes >= TIME_BIN(1) & relOnTimes < TIME_BIN(2);

    thisTbl.STIM_ON{t} = onTimes(inBin);
    thisTbl.conditions{t} = num2cell(pos(inBin,:), 2);
end

% Drop trials left with no qualifying flashes
thisTbl = thisTbl(cellfun(@numel, thisTbl.STIM_ON) > 0, :);
nStimuli = sum(cellfun(@numel, thisTbl.STIM_ON));
if isCellAlign
    alignStr = sprintf('%s(%d)', ALIGN_TO, ALIGN_IND);
else
    alignStr = ALIGN_TO;
end
fprintf('%d of %d correct trials have >=1 RF flash in [%d, %d] ms relative to %s (%d flashes total)\n', ...
    height(thisTbl), nTrialsTotal, TIME_BIN(1), TIME_BIN(2), alignStr, nStimuli);
if isempty(thisTbl)
    error('No RF flashes fall in [%d, %d] ms relative to %s -- nothing to plot.', TIME_BIN(1), TIME_BIN(2), alignStr);
end

end
