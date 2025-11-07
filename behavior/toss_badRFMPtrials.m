function [tbl_new] = toss_badRFMPtrials(tbl)
% TOSS_BADRFMPTRIALS Remove final RFMP stimulus for trials with broken fixation.
%
% This function takes a trial table `tbl` where each row is a single trial.
% The table must contain the columns:
%   - 'STIM_ON'   : cell array, each cell holding an array of stimulus-on times
%   - 'STIM_OFF'  : cell array, each cell holding an array of stimulus-off times
%   - 'conditions': cell array, each cell holding a list of condition labels
%   - 'result'    : categorical or char, trial outcome (e.g., 'CORRECT')
%
% For any trial where the monkey did not maintain fixation (i.e., result ~= 'CORRECT'),
% the final stimulus entry is removed from both STIM_ON and STIM_OFF. The conditions
% vector is then trimmed to match the new number of stimuli. Trials whose conditions
% become empty are removed entirely.

    % Remove final STIM_ON timestamp for non-correct trials
    tbl.STIM_ON(tbl.result ~= 'CORRECT') = ...
        cellfun(@(q) q(1:end-1), tbl.STIM_ON(tbl.result ~= 'CORRECT'), 'uni', 0);

    % Remove final STIM_OFF timestamp for non-correct trials
    tbl.STIM_OFF(tbl.result ~= 'CORRECT') = ...
        cellfun(@(q) q(1:end-1), tbl.STIM_OFF(tbl.result ~= 'CORRECT'), 'uni', 0);

    % Trim the conditions list so it matches the updated number of STIM_ON entries
    tbl.conditions = cellfun(@(conds, stimOn) conds(1:numel(stimOn)), ...
                             tbl.conditions, tbl.STIM_ON, 'uni', 0);

    % Drop any trials where conditions is now empty
    tbl_new = tbl(cellfun(@(q) ~isempty(q), tbl.conditions), :);
end