function p = paired_permutation_test(evoked_frs, baseline_frs, n_shuffles)
% Perform a paired permutation (sign-flip) test
%
% Inputs:
%   evoked_frs   - vector of firing rates during evoked epoch
%   baseline_frs - vector of firing rates during baseline epoch
%   n_shuffles     - number of permutations (e.g., 1000 or 5000)
%
% Output:
%   p              - two-sided p-value

    if nargin < 3
        n_shuffles = 1000;
    end

    % Ensure column vectors
    evoked_frs = evoked_frs(:);
    baseline_frs = baseline_frs(:);

    % Compute paired differences
    diff = evoked_frs - baseline_frs;

    % Observed statistic (mean difference)
    observed = mean(diff);

    % Preallocate null distribution
    null_dist = zeros(n_shuffles, 1);

    % Permutation via sign flipping
    for i = 1:n_shuffles
        signs = randi([0, 1], size(diff)) * 2 - 1;  % random +/-1
        null_dist(i) = mean(diff .* signs);
    end

    % Two-sided p-value
    p = mean(abs(null_dist) >= abs(observed));

end