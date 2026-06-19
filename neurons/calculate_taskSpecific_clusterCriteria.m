function [Tclust, Ttask] = calculate_taskSpecific_clusterCriteria(Tclust, Ttask, varargin)
%
% calculate_taskSpecific_clusterCriteria
%
% Computes epoch-specific quality and task metrics for a set of clusters,
% adding them as prefixed columns to Tclust. Call once per epoch (visual,
% saccade, pursuit) with the appropriate trial table and alignment parameters.
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
% -------------------------------------------------------------------------
% Tclust : cluster table, one row per unit (may span multiple probes).
%          Required columns:
%            cluster_id     - 0-based index into spiketimes cell arrays
%            probe_index    - which probe (matches spiketimes_<n> column)
%            probe_label    - hemisphere label ('lFEF' or 'rFEF')
%            unit_locations - [x y] position (µm)
%            probe_depth_mm - probe insertion depth (mm)
%
% Ttask  : trial table for the epoch of interest. Required columns:
%            spiketimes_<n>     - cell array of spike time vectors per cluster
%            END_TRIAL          - trial duration (ms)
%            angle              - target/motion direction (deg)
%            <RESPONSE_ALIGNTO> - alignment timestamps for response window
%            <BASELINE_ALIGNTO> - alignment timestamps for baseline window
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE PARAMETERS
% -------------------------------------------------------------------------
% 'EPOCH_NICKNAME'   char,  default 'sac'           Prefix for all output columns
% 'RESPONSE_WINDOW'  1×2,  default [-50 50]         Response epoch (ms)
% 'RESPONSE_ALIGNTO' char,  default 'saccadeOnset'  Alignment event for response
% 'BASELINE_WINDOW'  1×2,  default [-200 -100]      Baseline epoch (ms)
% 'BASELINE_ALIGNTO' char,  default 'saccadeOnset'  Alignment event for baseline
% 'CONDITIONS'       cell,  default {'angle','distance'}  Condition columns in Ttask
%
% -------------------------------------------------------------------------
% OUTPUT COLUMNS  (all prefixed by EPOCH_NICKNAME, e.g. 'sac_dp')
% -------------------------------------------------------------------------
% Presence / firing rate:
%   _trl_ratio       fraction of trials where unit fired at least once
%   _mean_fr         mean firing rate over full trial duration (Hz)
%
% Temporal stability  (epoch-agnostic, based on total per-trial FR):
%   _nBlocks         number of condition blocks detected
%   _blkPresRatio    min/max block presence ratio  (1 = active every trial)
%   _condSlope       mean fractional FR drift per repetition across conditions
%   _fracCP          fraction of conditions with a detected change point
%
% Epoch modulation:
%   _baseFR          mean baseline FR (Hz)
%   _mnFR            mean response FR (Hz)
%   _mnFR_peakWin    mean FR of the peak 25-ms sub-window within the response epoch (Hz)
%   _modIndex        response / baseline ratio
%   _pval_sr         Wilcoxon signed-rank p-value (raw, apply Bonferroni externally)
%   _exc             true if mean response > mean baseline (excitatory)
%   _dp              d′ (positive = excitatory, negative = suppressed)
%   _dpAbs           absolute d′
%   _dpPerCond       d′ per condition  {nClusters × 1 cell}
%   _pvPerCond       permutation p-value per condition
%   _pvalDP          true if any condition has d′ p < 0.05
%   _stdPerCond      response FR std per condition
%   _ffPerCond       Fano factor per condition
%
% Direction tuning:
%   _pvalDir         true if tuning is significant in any non-direction condition
%   _selDir          selectivity index  (0 = uniform, 1 = single direction)
%   _prefDir         vector-averaged preferred direction (deg)
%   _rhoLst          5th-percentile CI per direction
%   _rhoUst          95th-percentile CI per direction
%   _frs_perAng      per-trial FR distributions per direction  {cell}
%   _mnFR_perAng     mean FR per direction
%   _mxFR            peak mean FR across all directions
%   _bestDir         direction with highest mean FR (deg)
%   _hemi            'contra' or 'ipsi' relative to probe hemisphere

%% Input parsing
p = inputParser;
addRequired(p,  'Tclust');
addRequired(p,  'Ttask');
addParameter(p, 'EPOCH_NICKNAME',   'sac',                @ischar);
addParameter(p, 'RESPONSE_WINDOW',  [-50, 50],            @isnumeric);
addParameter(p, 'RESPONSE_ALIGNTO', 'saccadeOnset',       @ischar);
addParameter(p, 'BASELINE_WINDOW',  [-200, -100],         @isnumeric);
addParameter(p, 'BASELINE_ALIGNTO', 'saccadeOnset',       @ischar);
addParameter(p, 'CONDITIONS',       {'angle','distance'}, @iscell);
parse(p, Tclust, Ttask, varargin{:});

EPOCH_NICKNAME   = p.Results.EPOCH_NICKNAME;
RESPONSE_WINDOW  = p.Results.RESPONSE_WINDOW;
RESPONSE_ALIGNTO = p.Results.RESPONSE_ALIGNTO;
BASELINE_WINDOW  = p.Results.BASELINE_WINDOW;
BASELINE_ALIGNTO = p.Results.BASELINE_ALIGNTO;
CONDITIONS       = p.Results.CONDITIONS;

% Convert numeric alignment columns to cell arrays if needed
% (cellfun below uses u(1) to extract the scalar timestamp from each cell)
if ~iscell(Ttask.(RESPONSE_ALIGNTO)), Ttask.(RESPONSE_ALIGNTO) = num2cell(Ttask.(RESPONSE_ALIGNTO)); end
if ~iscell(Ttask.(BASELINE_ALIGNTO)), Ttask.(BASELINE_ALIGNTO) = num2cell(Ttask.(BASELINE_ALIGNTO)); end

probes = unique(Tclust.probe_index);

%% Remove trials where no unit fired (spike-sorting floor artifacts)
spikeCols = contains(Ttask.Properties.VariableNames, 'spiketimes');
n_empty   = sum(cellfun(@(v) sum(cellfun(@isempty, v)), Ttask{:, spikeCols}, 'uni', 1), 2);
Ttask     = Ttask(n_empty < height(Tclust), :);

%% Remove trials with population firing rate > 3 SD from session mean
% Compute response-window spike count matrix: [nTrials × nClusters]
pop_fr = nan(height(Ttask), height(Tclust));
for clust = 1:height(Tclust)
    spks = cellfun(@(q) q{Tclust.cluster_id(clust)+1}, ...
        Ttask.(sprintf('spiketimes_%d', Tclust.probe_index(clust))), 'uni', 0);
    pop_fr(:,clust) = cellfun(@(u,w) sum(w >= (u(1)+RESPONSE_WINDOW(1)) & ...
                                          w <  (u(1)+RESPONSE_WINDOW(2))), ...
        Ttask.(RESPONSE_ALIGNTO), spks, 'uni', 1);
end
pop_fr = pop_fr .* (1000 / diff(RESPONSE_WINDOW));

% Identify outlier trials: population mean FR > ±3 SD of the session distribution
trl_mean = mean(pop_fr, 2);      % mean across clusters per trial [nTrials × 1]
sess_mu  = mean(trl_mean);
sess_sd  = std(trl_mean);
Ttask    = Ttask(abs(trl_mean - sess_mu) <= 3*sess_sd, :);

%% Cluster depth relative to brain surface  (negative = deeper)
Tclust.cluster_depth_mm = ...
    (cellfun(@(q) q(2), Tclust.unit_locations, 'uni', 1) - Tclust.probe_depth_mm*1000) ./ 1000;

%% Presence ratio, mean FR, and temporal stability
% All three require per-cluster iteration, so they share a single loop.
%
% Stability uses total per-trial FR (epoch-agnostic) and the block structure
% derived from the condition sequence (a new block begins when any condition
% repeats). No criteria are applied here — raw values stored for later use.

blk_ids   = assign_blocks(Ttask{:, CONDITIONS});
blk_conds = Ttask{:, CONDITIONS};

[Tclust.([EPOCH_NICKNAME, '_nBlocks']),      ...
 Tclust.([EPOCH_NICKNAME, '_blkPresRatio']), ...
 Tclust.([EPOCH_NICKNAME, '_condSlope']),    ...
 Tclust.([EPOCH_NICKNAME, '_fracCP'])]       = deal(nan(height(Tclust), 1));

for clust = 1:height(Tclust)
    spk_col = sprintf('spiketimes_%d', Tclust.probe_index(clust));
    col     = Tclust.cluster_id(clust) + 1;
    spks    = cellfun(@(q) q{col}, Ttask.(spk_col), 'uni', 0);
    fr_trl  = cellfun(@(q,t) numel(q)/t*1000, spks, num2cell(Ttask.END_TRIAL), 'uni', 1);

    Tclust.([EPOCH_NICKNAME, '_trl_ratio'])(clust) = mean(cellfun(@(q) ~isempty(q), spks));
    Tclust.([EPOCH_NICKNAME, '_mean_fr'])(clust)   = mean(fr_trl);

    [Tclust.([EPOCH_NICKNAME, '_blkPresRatio'])(clust), ...
     Tclust.([EPOCH_NICKNAME, '_nBlocks'])(clust)]     = block_presence_ratio(fr_trl, blk_ids);
    [Tclust.([EPOCH_NICKNAME, '_condSlope'])(clust),   ...
     Tclust.([EPOCH_NICKNAME, '_fracCP'])(clust)]       = condition_stability(fr_trl, blk_conds);
end

%% Baseline and response firing rates  [nTrials × nClusters, Hz]
[base_fr, resp_fr] = deal(nan(height(Ttask), height(Tclust)));
for clust = 1:height(Tclust)
    spk_col = sprintf('spiketimes_%d', Tclust.probe_index(clust));
    col     = Tclust.cluster_id(clust) + 1;
    spks    = cellfun(@(q) q{col}, Ttask.(spk_col), 'uni', 0);

    base_fr(:,clust) = cellfun(@(u,w) sum(w >= (u(1)+BASELINE_WINDOW(1)) & w < (u(1)+BASELINE_WINDOW(2))), ...
        Ttask.(BASELINE_ALIGNTO), spks, 'uni', 1);
    resp_fr(:,clust) = cellfun(@(u,w) sum(w >= (u(1)+RESPONSE_WINDOW(1)) & w < (u(1)+RESPONSE_WINDOW(2))), ...
        Ttask.(RESPONSE_ALIGNTO), spks, 'uni', 1);
end
base_fr = base_fr .* (1000 / abs(diff(BASELINE_WINDOW)));
resp_fr = resp_fr .* (1000 / diff(RESPONSE_WINDOW));

%% Wilcoxon signed-rank test: response vs baseline
% Stores raw (uncorrected) p-values. Apply Bonferroni/Holm across the three
% epochs (visual, saccade, pursuit) externally after all epochs are computed.
pval_sr = nan(height(Tclust), 1);
for clust = 1:height(Tclust)
    [pval_sr(clust), ~] = signrank(resp_fr(:,clust), base_fr(:,clust));
end
Tclust.([EPOCH_NICKNAME, '_pval_sr']) = pval_sr;
Tclust.([EPOCH_NICKNAME, '_exc'])     = (mean(resp_fr)' > mean(base_fr)');

%% Peak 25-ms sub-window firing rate  [nClusters × 1, Hz]
% Slide a 25-ms window across the response epoch in the trial-averaged
% response, then take the maximum. Converted to Hz.
PEAK_WIN   = 25;  % ms
win_starts = RESPONSE_WINDOW(1) : (RESPONSE_WINDOW(2) - PEAK_WIN);
peak_fr    = nan(height(Tclust), 1);
for clust = 1:height(Tclust)
    spk_col = sprintf('spiketimes_%d', Tclust.probe_index(clust));
    col     = Tclust.cluster_id(clust) + 1;
    spks    = cellfun(@(q) q{col}, Ttask.(spk_col), 'uni', 0);

    mean_count_per_win = arrayfun(@(s) mean(cellfun(@(u, w) ...
        sum(w >= (u(1)+s) & w < (u(1)+s+PEAK_WIN)), ...
        Ttask.(RESPONSE_ALIGNTO), spks)), win_starts);
    peak_fr(clust) = max(mean_count_per_win) * (1000 / PEAK_WIN);
end

%% Epoch summary FRs and modulation index
% Replace exact zeros before division (modIndex only; does not affect d′ or signrank)
base_fr(base_fr == 0) = 1e-6;
resp_fr(resp_fr == 0) = 1e-6;

Tclust.([EPOCH_NICKNAME, '_baseFR'])      = mean(base_fr)';
Tclust.([EPOCH_NICKNAME, '_mnFR'])        = mean(resp_fr)';
Tclust.([EPOCH_NICKNAME, '_mnFR_peakWin']) = peak_fr;
Tclust.([EPOCH_NICKNAME, '_modIndex'])    = mean(resp_fr)' ./ mean(base_fr)';

%% Discriminability (d′) between response and baseline

% Overall d′ across all conditions and trials
dp = (mean(resp_fr,1) - mean(base_fr,1)) ./ sqrt(0.5*(var(resp_fr,0) + var(base_fr,0)));
Tclust.([EPOCH_NICKNAME, '_dp'])    = dp';
Tclust.([EPOCH_NICKNAME, '_dpAbs']) = (mean(abs(resp_fr-base_fr),1) ./ ...
    sqrt(0.5*(var(resp_fr,0) + var(base_fr,0))))';

% Per-condition d′, permutation p-value, variability, and Fano factor
condVals     = Ttask{:, CONDITIONS};
unique_conds = unique(condVals, 'rows');
n_conds      = size(unique_conds, 1);

[dp_perCond, pv_perCond, std_perCond, ff_perCond] = deal(nan(height(Tclust), n_conds));
for cond = 1:n_conds
    mask = ismember(condVals, unique_conds(cond,:), 'rows');
    bfr  = base_fr(mask, :);
    rfr  = resp_fr(mask, :);

    [dp_obs, p_val]     = cellfun(@(v,b) compute_dprime_perm(v,b), num2cell(rfr,1), num2cell(bfr,1), 'uni', 1);
    dp_perCond(:,cond)  = dp_obs';
    pv_perCond(:,cond)  = p_val';
    std_perCond(:,cond) = std(rfr)';
    ff_perCond(:,cond)  = (var(rfr) ./ mean(rfr))';
end

Tclust.([EPOCH_NICKNAME, '_dpPerCond'])  = num2cell(dp_perCond,  2);
Tclust.([EPOCH_NICKNAME, '_pvPerCond'])  = num2cell(pv_perCond,  2);
Tclust.([EPOCH_NICKNAME, '_pvalDP'])     = min(pv_perCond, [], 2) < 0.05;
Tclust.([EPOCH_NICKNAME, '_stdPerCond']) = num2cell(std_perCond, 2);
Tclust.([EPOCH_NICKNAME, '_ffPerCond'])  = num2cell(ff_perCond,  2);

%% Direction tuning
% Tuning significance: tested separately per non-direction condition (e.g., per amplitude
% or speed level); a unit passes if ANY level shows significant tuning.
% Selectivity and preferred direction are computed across all trials pooled.

nonDir_conds = CONDITIONS(~ismember(CONDITIONS, 'angle'));
condVals     = Ttask{:, nonDir_conds};
unique_conds = unique(condVals, 'rows');

pval_perCond = nan(height(Tclust), numel(unique_conds));
for cond = 1:numel(unique_conds)
    pvals = cell(numel(probes), 1);
    for prb = 1:numel(probes)
        [~,~,~,~,~, pval_dir] = calculate_tuning( ...
            Ttask(ismember(condVals, unique_conds(cond,:), 'rows'), :), ...
            'PROBE_INDEX', probes(prb), 'FR_WIN', RESPONSE_WINDOW, 'ALIGN_TO', RESPONSE_ALIGNTO);
        pvals{prb} = pval_dir;
    end
    pval_perCond(:,cond) = vertcat(pvals{:});
end
Tclust.([EPOCH_NICKNAME, '_pvalDir']) = any(pval_perCond, 2);

% Tuning curves, preferred direction, and hemisphere assignment
uniqueAngs = sort(unique(Ttask.angle));
[sel_dir, pref_dir, rhoLst, rhoUst, frs_per_ang, hemi] = deal(cell(numel(probes), 1));

for prb = 1:numel(probes)
    [sld, pd, rl, ru, fpa, ~] = calculate_tuning(Ttask, 'WITH_MAXSTAT', false, ...
        'PROBE_INDEX', probes(prb), 'FR_WIN', RESPONSE_WINDOW, 'ALIGN_TO', RESPONSE_ALIGNTO);

    % Hemisphere label for this probe (contra vs. ipsi based on preferred direction)
    prb_label = Tclust.probe_label(find(Tclust.probe_index == probes(prb), 1));
    is_right  = (pd >= 0 & pd < 90) | (pd > 270);

    sel_dir{prb}     = sld;
    pref_dir{prb}    = pd;
    rhoLst{prb}      = rl;
    rhoUst{prb}      = ru;
    frs_per_ang{prb} = fpa';
    hemi{prb}        = assign_hemi(is_right, prb_label);
end

Tclust.([EPOCH_NICKNAME, '_selDir'])     = vertcat(sel_dir{:});
Tclust.([EPOCH_NICKNAME, '_prefDir'])    = vertcat(pref_dir{:});
Tclust.([EPOCH_NICKNAME, '_rhoLst'])     = vertcat(rhoLst{:});
Tclust.([EPOCH_NICKNAME, '_rhoUst'])     = vertcat(rhoUst{:});
Tclust.([EPOCH_NICKNAME, '_frs_perAng']) = vertcat(frs_per_ang{:});
Tclust.([EPOCH_NICKNAME, '_hemi'])       = vertcat(hemi{:});

% Per-direction summary stats
Tclust.([EPOCH_NICKNAME, '_mnFR_perAng']) = cellfun(@(q) mean(q), ...
    Tclust.([EPOCH_NICKNAME, '_frs_perAng']), 'uni', 1);

[mxFR, mxIdx] = max(Tclust.([EPOCH_NICKNAME, '_mnFR_perAng']), [], 2);
Tclust.([EPOCH_NICKNAME, '_mxFR'])    = mxFR;
Tclust.([EPOCH_NICKNAME, '_bestDir']) = uniqueAngs(mxIdx);

end


% =========================================================================
function block_ids = assign_blocks(conditions)
% Assign trials to sequential blocks. A new block begins whenever any
% condition repeats within the current block.
    n         = size(conditions, 1);
    block_ids = ones(n, 1);
    block_num = 1;
    seen      = conditions(1, :);
    for t = 2:n
        cond = conditions(t, :);
        if any(all(seen == cond, 2))
            block_num = block_num + 1;
            seen      = cond;
        else
            seen      = [seen; cond];
        end
        block_ids(t) = block_num;
    end
end


% =========================================================================
function [ratio, n_blocks] = block_presence_ratio(fr_series, block_ids)
% Returns the min/max block presence ratio: the fraction of trials in
% which the unit fired at least once, taken as minimum across blocks.
% A value of 1 means the unit fired on every trial in every block.
    blocks   = unique(block_ids);
    n_blocks = numel(blocks);
    pr       = nan(n_blocks, 1);
    for b = 1:n_blocks
        pr(b) = mean(fr_series(block_ids == blocks(b)) > 0);
    end
    ratio = min(pr, [], 'omitnan') / max(max(pr, [], 'omitnan'), eps);
end


% =========================================================================
function [mean_slope, frac_cp] = condition_stability(fr_series, conditions)
% Assesses firing rate stability across repetitions of each condition.
%
% mean_slope : mean fractional FR change per repetition (0 = flat)
% frac_cp    : fraction of conditions containing a detected change point
%
% Constants below were chosen to minimize false positives on short
% per-condition sequences; tune after inspecting population distributions.
    SMOOTH_WIN = 3;    % repetitions to smooth before change-point detection
    MIN_REPS   = 4;    % skip conditions with fewer repetitions than this
    CP_THRESH  = 5;    % minimum FR shift (Hz) to count as a change point

    unique_conds = unique(conditions, 'rows');
    n_conds      = size(unique_conds, 1);
    slopes       = nan(n_conds, 1);
    has_cp       = false(n_conds, 1);

    for c = 1:n_conds
        cond_mask = all(conditions == unique_conds(c,:), 2);
        fr_cond   = fr_series(cond_mask);

        if numel(fr_cond) < MIN_REPS, continue; end

        % Normalized slope: fractional FR change per repetition
        mn_fr = mean(fr_cond, 'omitnan');
        if mn_fr > 0
            x         = (1:numel(fr_cond))';
            b         = [ones(size(x)), x] \ fr_cond;
            slopes(c) = b(2) / mn_fr;
        end

        % Change-point detection on smoothed series
        % MinThreshold derived from CP_THRESH (Hz) and sequence length
        if numel(fr_cond) >= SMOOTH_WIN * 2
            smoothed   = movmean(fr_cond, SMOOTH_WIN);
            min_thresh = (CP_THRESH^2) * numel(smoothed) / 4;
            ipt        = findchangepts(smoothed', 'Statistic', 'mean', ...
                                       'MinThreshold', min_thresh);
            has_cp(c)  = ~isempty(ipt);
        end
    end

    mean_slope = mean(slopes, 'omitnan');
    frac_cp    = mean(has_cp, 'omitnan');
end


% =========================================================================
function hemi = assign_hemi(is_right_dir, probe_label)
% Map preferred direction laterality to contra/ipsi given probe hemisphere.
    hemi = strings(numel(is_right_dir), 1);
    if probe_label == "lFEF"
        hemi(is_right_dir)  = "contra";
        hemi(~is_right_dir) = "ipsi";
    else
        hemi(is_right_dir)  = "ipsi";
        hemi(~is_right_dir) = "contra";
    end
    hemi = categorical(hemi);
end
