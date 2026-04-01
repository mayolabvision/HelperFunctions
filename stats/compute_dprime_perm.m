function [dp_obs, p_val] = compute_dprime_perm(v, b, n_shuffles)
% Computes d' and permutation-based p-value for one neuron
%
% Inputs:
%   v       - response firing rates (trials x 1)
%   b       - baseline firing rates (trials x 1)
%   n_perm  - number of permutations (e.g., 1000)
%
% Outputs:
%   dp_obs  - observed d-prime
%   p_val   - two-sided p-value from permutation test
%   dp_null - null distribution of d-prime (optional, for debugging/plots)

    if nargin < 3
        n_shuffles = 1000;
    end

    % Remove NaNs just in case
    v = v(~isnan(v));
    b = b(~isnan(b));

    % Observed d'
    denom = sqrt(0.5 * (var(v,0) + var(b,0)));
    if denom == 0
        dp_obs = NaN;
        p_val = NaN;
        return
    end

    dp_obs = (mean(v) - mean(b)) / denom;

    % Permutation setup
    all_vals = [v; b];
    n_v = length(v);
    n_b = length(b);

    dp_null = zeros(n_shuffles,1);

    for i = 1:n_shuffles
        idx = randperm(n_v + n_b);

        v_shuff = all_vals(idx(1:n_v));
        b_shuff = all_vals(idx(n_v+1:end));

        denom_shuff = sqrt(0.5 * (var(v_shuff,0) + var(b_shuff,0)));

        if denom_shuff == 0
            dp_null(i) = NaN;
        else
            dp_null(i) = (mean(v_shuff) - mean(b_shuff)) / denom_shuff;
        end
    end

    % Remove NaNs from null
    dp_null = dp_null(~isnan(dp_null));

    % Two-sided p-value
    p_val = mean(abs(dp_null) >= abs(dp_obs));

end