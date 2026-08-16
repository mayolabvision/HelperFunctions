function [mm_frs, mm_indices] = meanMatch_pairsRSC(frCells, varargin)
% meanMatch_pairsRSC
%
% Perform firing-rate mean matching across conditions by iteratively
% removing neuron pairs until the mean firing rate within each condition
% is within a specified percentage of the population mean firing rate.
% Pair removal can be deterministic, random, or biased-random, and is
% applied identically to firing-rate and spike-count correlation data.
%
% INPUTS
% ------
% frCells : cell array
%     Cell array of size [N x M] (or any shape), where each cell contains
%     a numeric vector of mean firing rates for neuron pairs belonging to
%     a specific condition (e.g., d' × RT decile). Each element of the
%     vector corresponds to one neuron pair.
%
% rscCells : cell array
%     Cell array of the same size as frCells. Each cell contains a numeric
%     vector of spike-count correlation values (r_sc) for the same neuron
%     pairs and in the same order as the corresponding frCells entry.
%
% OPTIONAL PARAMETERS (Name–Value pairs)
% --------------------------------------
% 'PERCENT_WITHIN' : scalar (default = 5)
%     Percent tolerance for mean matching. Matching succeeds when the
%     mean firing rate within a cell lies within ±PERCENT_WITHIN percent
%     of the population mean firing rate.
%
% 'MIN_PAIRS' : scalar (default = 5)
%     Minimum number of neuron pairs that must remain in a cell after
%     pair removal. If this threshold is reached before matching succeeds,
%     the cell is returned empty.
%
% 'DROP_MODE' : string (default = 'directed')
%     Strategy used to remove neuron pairs:
%       - 'directed'      : deterministically removes the largest firing
%                           rate if the mean is too high, or the smallest
%                           firing rate if the mean is too low.
%       - 'random'        : removes a randomly selected pair.
%       - 'biased-random' : removes a randomly selected pair from the
%                           upper or lower tail of the firing-rate
%                           distribution, depending on whether the mean
%                           is above or below the target range.
%
% 'BIAS_FRAC' : scalar in (0,1) (default = 0.10)
%     Fraction of pairs defining the upper or lower tail used for
%     biased-random removal. For example, BIAS_FRAC = 0.10 selects from
%     the top or bottom 10% of firing rates.
%
% OUTPUTS
% -------
% mm_frs : cell array
%     Cell array of the same size as frCells, containing mean-matched
%     firing-rate vectors after pair removal. Cells for which matching
%     failed are returned empty.
%
% mm_rsc : cell array
%     Cell array of the same size as rscCells, containing spike-count
%     correlation vectors corresponding to the retained neuron pairs.
%
% NOTES
% -----
% - Mean matching is performed independently for each cell.
% - Pair removal is applied identically to firing rates and r_sc values.
% - This function does not perform repeated resampling; stochastic modes
%   ('random' or 'biased-random') should be repeated externally if
%   averaging across repetitions is desired.

p = inputParser;
addParameter(p, 'PERCENT_WITHIN', 1, @isnumeric);
addParameter(p, 'MIN_PAIRS', 100, @isnumeric);
addParameter(p, 'DROP_MODE', 'directed', ...
    @(x) any(validatestring(x, {'random','directed','biased-random'})));
addParameter(p, 'BIAS_FRAC', 0.25, @isnumeric);

parse(p,varargin{:});

PERCENT_WITHIN  =  p.Results.PERCENT_WITHIN;
MIN_PAIRS       =  p.Results.MIN_PAIRS;
DROP_MODE       =  p.Results.DROP_MODE;
BIAS_FRAC       =  p.Results.BIAS_FRAC;

%=======================================================================

% Population mean firing rate
allFRs = vertcat(frCells{:});
pop_mnFR = mean(allFRs);

mnFR_targetRange = pop_mnFR * ...
    [1 - PERCENT_WITHIN/100, 1 + PERCENT_WITHIN/100];

%%=======================================================================

[mm_frs, mm_indices] = cellfun(@(fr) mean_match(fr, mnFR_targetRange, ...
                                     MIN_PAIRS, DROP_MODE, BIAS_FRAC), ...
                                     frCells, 'uni', 0);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mm_frs, mm_indices] = mean_match(mm_frs, mnFR_targetRange, ...
                                           MIN_PAIRS, DROP_MODE, BIAS_FRAC)

totalPairs = numel(mm_frs);
dropped_total = 0;

original_frs = mm_frs;       % store original values
original_indices = 1:totalPairs;

switch DROP_MODE
    case 'random'
        success = false;
        while ~success
            mm_frs = original_frs;
            mm_indices = original_indices;
            dropped = 0;

            while numel(mm_frs) > MIN_PAIRS
                mnFR0 = mean(mm_frs);

                % Check if we're within target
                if mnFR0 >= mnFR_targetRange(1) && mnFR0 <= mnFR_targetRange(2)
                    success = true;
                    dropped_total = dropped_total + dropped;
                    fprintf('%d/%d pairs dropped\n', dropped_total, totalPairs);
                    return
                end

                % Drop a random pair
                idx = randi(numel(mm_frs));
                mm_frs(idx) = [];
                mm_indices(idx) = [];
                dropped = dropped + 1;
            end

            % If we reach MIN_PAIRS and still not in range, restart
            fprintf('Random attempt failed; restarting...\n');
        end

    otherwise
        % Non-random modes keep original logic
        mm_indices = 1:totalPairs;
        dropped = 0;

        while numel(mm_frs) > MIN_PAIRS
            mnFR0 = mean(mm_frs);

            % Check if we're within target
            if mnFR0 >= mnFR_targetRange(1) && mnFR0 <= mnFR_targetRange(2)
                return
            end

            % Decide which pair to drop
            switch DROP_MODE
                case 'directed'
                    if mnFR0 > mnFR_targetRange(2)
                        [~, idx] = max(mm_frs);
                    else
                        [~, idx] = min(mm_frs);
                    end

                case 'biased-random'
                    n = numel(mm_frs);
                    k = max(1, ceil(BIAS_FRAC * n)); % top/bottom percent

                    if mnFR0 > mnFR_targetRange(2)
                        % mean too high → drop from upper tail
                        [~, ord] = sort(mm_frs, 'descend');
                    else
                        % mean too low → drop from lower tail
                        [~, ord] = sort(mm_frs, 'ascend');
                    end

                    candidates = ord(1:k);
                    idx = candidates(randi(numel(candidates)));
            end

            % Drop selected pair
            mm_frs(idx) = [];
            mm_indices(idx) = [];
            dropped = dropped + 1;
        end

        % If we get here, matching failed
        fprintf('Mean match failed (%d/%d pairs kept)\n', numel(mm_frs), totalPairs);
        mm_frs = [];
        mm_indices = [];
end

end