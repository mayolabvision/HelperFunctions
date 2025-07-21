function [np_mask, ripple_mask] = match_syncPulses_RipToNP(np_pulse_timeStamps, ripple_pulse_timeStamps)
% match_syncPulses_RipToNP - Match pulse timestamps between NP and Ripple.
%
% Inputs:
%   np_pulse_timeStamps     - Nx1 datetime array (all NP pulses)
%   ripple_pulse_timeStamps - Mx1 datetime array (subset of "good" pulses from Ripple)
%
% Outputs:
%   np_mask                 - Nx1 logical array, true where a match to Ripple is found
%   ripple_mask             - Mx1 logical array, true where a match to NP is found

    tol_ms = 20;
    tol = milliseconds(tol_ms);

    np_mask = false(size(np_pulse_timeStamps));
    ripple_mask = false(size(ripple_pulse_timeStamps));

    for i = 1:numel(ripple_pulse_timeStamps)
        time_diffs = abs(np_pulse_timeStamps - ripple_pulse_timeStamps(i));
        [min_diff, min_idx] = min(time_diffs);
        if min_diff <= tol
            np_mask(min_idx) = true;
            ripple_mask(i) = true;
        end
    end

end