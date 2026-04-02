function [Tclust, Ttask] = calculate_taskSpecific_clusterCriteria(Tclust, Ttask, varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
% -------------------------------------------------------------------------
% Tclust : table
%     One row per cluster (unit). May contain clusters recorded across
%     multiple probes. Must include at least:
%         - cluster_id        : index of cluster within probe
%         - probe_index       : probe number corresponding to spiketimes_X
%         - unit_locations    : [x y] position of unit (µm)
%         - probe_depth_mm    : insertion depth of probe (mm)
%
% Ttask : table
%     One row per trial. Should already exclude trials that should not be
%     analyzed (e.g., incorrect trials or behavioral exclusions).
%
%     Required columns:
%         - spiketimes_<probe> : cell array of spike times per cluster
%         - END_TRIAL           : trial duration (ms)
%         - alignment variables specified below
%
% -------------------------------------------------------------------------
% OPTIONAL NAME–VALUE PARAMETERS
% -------------------------------------------------------------------------
%
% 'EPOCH_NICKNAME' (char, default = 'sac')
%     Prefix used when adding computed metrics to Tclust, can be whatever.
%     Example output variables:
%         sac_mean_fr, sac_dp, sac_prefDir, etc.
%
% 'RESPONSE_WINDOW' (1x2 double, default = [-50 50])
%     Time window (ms) used to compute response firing rate.
%
% 'RESPONSE_ALIGNTO' (char, default = 'saccadeOnset')
%     Column in Ttask containing alignment timestamps for the response window.
%
% 'BASELINE_WINDOW' (1x2 double, default = [-200 -100])
%     Time window (ms) used to compute baseline firing rate.
%
% 'BASELINE_ALIGNTO' (char, default = 'saccadeOnset')
%     Column in Ttask containing alignment timestamps for the baseline window.
%
% 'CONDITIONS' (cell array of char, default = {'angle','distance'})
%     Column names in Ttask defining experimental conditions. Unique
%     combinations are used for per-condition statistics.
%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------
% 1. Removes trials where no neurons fired (sorting cutoff artifacts).
% 2. Rejects trials whose population firing rate is >3 SD from mean.
% 3. Computes per-cluster:
%       - depth relative to brain surface
%       - trial participation ratio
%       - overall mean firing rate
% 4. Computes baseline and response firing rates per trial.
% 5. Computes discriminability metrics:
%       - d′ between baseline and response
%       - absolute d′
%       - permutation-test d′ per condition
%       - Fano factor and variability per condition
% 6. Computes direction tuning metrics:
%       - tuning significance
%       - selectivity index
%       - preferred direction
%       - confidence intervals
%       - firing rate per direction
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% Tclust : table
%     Input cluster table augmented with task-specific metrics whose names
%     are prefixed by EPOCH_NICKNAME.
%
% Ttask : table
%     Trial table after removal of excluded trials.
%

p = inputParser;
addRequired(p, 'Tclust', @istable);
addRequired(p, 'Ttask', @istable);

addParameter(p, 'EPOCH_NICKNAME', 'sac', @ischar);
addParameter(p, 'RESPONSE_WINDOW', [-50,50], @isnumeric);
addParameter(p, 'RESPONSE_ALIGNTO', 'saccadeOnset', @ischar);
addParameter(p, 'BASELINE_WINDOW', [-200,-100], @isnumeric);
addParameter(p, 'BASELINE_ALIGNTO', 'saccadeOnset', @ischar);
addParameter(p, 'CONDITIONS', {'angle','distance'}, @iscell);

parse(p, Tclust, Ttask, varargin{:});

EPOCH_NICKNAME    =  p.Results.EPOCH_NICKNAME;
RESPONSE_WINDOW   =  p.Results.RESPONSE_WINDOW;
RESPONSE_ALIGNTO  =  p.Results.RESPONSE_ALIGNTO;
BASELINE_WINDOW   =  p.Results.BASELINE_WINDOW;
BASELINE_ALIGNTO  =  p.Results.BASELINE_ALIGNTO;
CONDITIONS        =  p.Results.CONDITIONS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(Ttask.(RESPONSE_ALIGNTO)), Ttask.(RESPONSE_ALIGNTO) = num2cell(Ttask.(RESPONSE_ALIGNTO)); end
if ~iscell(Ttask.(BASELINE_ALIGNTO)), Ttask.(BASELINE_ALIGNTO) = num2cell(Ttask.(BASELINE_ALIGNTO)); end

% Number of unique probes, for example if you recorded dual-hemi
probes = unique(Tclust.probe_index);

%% Removes trials where no neurons fired (sorting cutoff artifacts)

spikeCols = contains(Ttask.Properties.VariableNames, 'spiketimes');
Ttask = Ttask(sum(cellfun(@(v) sum(cellfun(@(q) isempty(q), v, 'uni', 1)), Ttask{:,spikeCols}, 'uni', 1) == height(Ttask),2)==0,:);

%% Rejects trials whose population firing rate is >3 SD from mean.

% resp_fr --> trials x clusters, in Hz
resp_fr = nan(height(Ttask),height(Tclust));
for clust = 1:height(Tclust)
    spks = cellfun(@(q) q{Tclust.cluster_id(clust)+1}, Ttask.(sprintf('spiketimes_%d',Tclust.probe_index(clust))), 'uni', 0);
    resp_fr(:,clust) = cellfun(@(u,w) sum(w>=(u(1)+RESPONSE_WINDOW(1)) & w<(u(1)+RESPONSE_WINDOW(2))), Ttask.(RESPONSE_ALIGNTO), spks, 'uni', 1);
end
resp_fr = resp_fr .* (1000/(RESPONSE_WINDOW(2)-RESPONSE_WINDOW(1)));

resp_mean = mean(resp_fr,2);
resp_std = std(resp_fr,[],2);
resp_idx = resp_mean < (resp_mean - resp_std .* 3) | resp_mean > (resp_mean + resp_std .* 3);

Ttask = Ttask(~resp_idx, :);

%% ==== Calculates depth relative to brain surface, trial/presence ratio, and mean firing rate ====

% cluster_depth_mm  =  depth of cluster, relative to brain surface.... more negative = deeper 
Tclust.cluster_depth_mm = (cellfun(@(q) q(2), Tclust.unit_locations, 'uni', 1) - (Tclust.probe_depth_mm*1000))./1000;

% trl_ratio  =  proportion of trials cluster fired at least once
% mean_fr    =  firing rate of cluster within full duration of trials, in Hz
for clust = 1:height(Tclust) 
    Tclust.([EPOCH_NICKNAME, '_trl_ratio'])(clust) = sum(cellfun(@(q) ~isempty(q), [cellfun(@(q) q{Tclust.cluster_id(clust)+1}, Ttask.(sprintf('spiketimes_%d',Tclust.probe_index(clust))), 'uni', 0)], 'uni', 1)) / (height(Ttask)); 
    Tclust.([EPOCH_NICKNAME, '_mean_fr'])(clust) = mean(cellfun(@(q,t) ((numel(q{Tclust.cluster_id(clust)+1}))./t)*1000, Ttask.(sprintf('spiketimes_%d',Tclust.probe_index(clust))), num2cell(Ttask.END_TRIAL), 'uni', 1));
end

%% ==== Computes baseline and response firing rates per cluster and trial ===

% base_fr, resp_fr --> trials x clusters, in Hz
[base_fr, resp_fr] = deal(nan(height(Ttask),height(Tclust)));
for clust = 1:height(Tclust)
    spks = cellfun(@(q) q{Tclust.cluster_id(clust)+1}, Ttask.(sprintf('spiketimes_%d',Tclust.probe_index(clust))), 'uni', 0);

    base_fr(:,clust) = cellfun(@(u,w) sum(w>=(u(1)+BASELINE_WINDOW(1)) & w<(u(1)+BASELINE_WINDOW(2))), Ttask.(BASELINE_ALIGNTO), spks, 'uni', 1);
    resp_fr(:,clust) = cellfun(@(u,w) sum(w>=(u(1)+RESPONSE_WINDOW(1)) & w<(u(1)+RESPONSE_WINDOW(2))), Ttask.(RESPONSE_ALIGNTO), spks, 'uni', 1);
end
base_fr = base_fr .* (1000/(BASELINE_WINDOW(2)-BASELINE_WINDOW(1)));
resp_fr = resp_fr .* (1000/(RESPONSE_WINDOW(2)-RESPONSE_WINDOW(1)));

base_fr(base_fr==0) = 1e-6; resp_fr(resp_fr==0) = 1e-6;

% baseFR = mean firing rate during baseline epoch, across all conditions
Tclust.([EPOCH_NICKNAME, '_baseFR']) = mean(base_fr)';

% mnFR   = mean firing rate during response epoch, across all conditions
Tclust.([EPOCH_NICKNAME, '_mnFR']) = mean(resp_fr)';

% modIndex = modulation index between response / baseline
Tclust.([EPOCH_NICKNAME, '_modIndex']) = mean(resp_fr)' ./ mean(base_fr)';

%% ==== Computes discriminability metrics between response and baseline windows ====

% dp = discriminability between baseline and response activity, positive = FR tends to increase and negative = FR tends to decrease
dp = (mean(resp_fr,1) - mean(base_fr,1)) ./ (sqrt(0.5*(var(resp_fr,0) + var(base_fr,0))));
Tclust.([EPOCH_NICKNAME, '_dp']) = dp';

% dpAbs = discriminability between baseline and response activity, ignoring if FR increases/decreases on different trials
dp_abs = mean(abs(resp_fr-base_fr),1) ./ sqrt(0.5*(var(resp_fr,0) + var(base_fr,0)));
Tclust.([EPOCH_NICKNAME, '_dpAbs']) = dp_abs';

% dp and fano factor per unique condition
condVals = Ttask{:, CONDITIONS};
unique_conds = unique(condVals,'rows');
[dp_perCond,pv_perCond,std_perCond,ff_perCond] = deal(nan(height(Tclust),size(unique_conds,1)));
for cond = 1:size(unique_conds,1)
    bfr = base_fr(ismember(condVals, unique_conds(cond,:), 'rows'),:);
    rfr = resp_fr(ismember(condVals, unique_conds(cond,:), 'rows'),:);

    [dp_obs, p_val] = cellfun(@(v,b) compute_dprime_perm(v,b), num2cell(rfr,1), num2cell(bfr,1), 'uni', 1);
    dp_perCond(:,cond) = dp_obs';
    pv_perCond(:,cond) = p_val';

    std_perCond(:,cond) = std(rfr)';
    ff_perCond(:,cond) = var(rfr)' ./ mean(rfr)';
end

% dpPerCond = d-prime per unique condition
Tclust.([EPOCH_NICKNAME, '_dpPerCond']) = num2cell(dp_perCond,2);

% pvPerCond = p-value of d-prime per unique condition
Tclust.([EPOCH_NICKNAME, '_pvPerCond']) = num2cell(pv_perCond,2);

% pvalDP = logical index if at least one condition has significant p-value
Tclust.([EPOCH_NICKNAME, '_pvalDP']) = min(pv_perCond,[],2)<0.05;

% stdPerCond = FR std per unique condition
Tclust.([EPOCH_NICKNAME, '_stdPerCond']) = num2cell(std_perCond,2);

% stdPerCond = fano factor per unique condition
Tclust.([EPOCH_NICKNAME, '_ffPerCond']) = num2cell(ff_perCond,2);

%% ==== Computes direction tuning metrics ====

nonDir_conds = CONDITIONS(~ismember(CONDITIONS,'angle'));
condVals = Ttask{:, nonDir_conds};
unique_conds = unique(condVals,'rows');

pval_perCond = nan(height(Tclust),numel(unique_conds));
for cond = 1:numel(unique_conds)
    pvals = cell(numel(probes),1);
    for prb = 1:numel(probes)
        [~, ~, ~, ~, ~, pval_dir] = calculate_tuning(Ttask(ismember(condVals, unique_conds(cond,:), 'rows'),:), ...
                                                     'PROBE_INDEX', probes(prb), 'FR_WIN', RESPONSE_WINDOW, 'ALIGN_TO', RESPONSE_ALIGNTO);
        pvals{prb} = pval_dir;
    end
    pval_perCond(:,cond) = vertcat(pvals{:});
end

% pvalDir = logical value for if any of the conditions (amplitude or speed) have significant direction tuning
Tclust.([EPOCH_NICKNAME, '_pvalDir']) = any(pval_perCond, 2);

[sel_dir, pref_dir, rhoLst, rhoUst, frs_per_ang] = deal(cell(numel(probes),1));
for prb = 1:numel(probes)
    [sld, pd, rl, ru, fpa, ~] = calculate_tuning(Ttask, 'WITH_MAXSTAT', false, ...
                                                 'PROBE_INDEX', probes(prb), 'FR_WIN', RESPONSE_WINDOW, 'ALIGN_TO', RESPONSE_ALIGNTO);

    sel_dir{prb} = sld; pref_dir{prb} = pd; rhoLst{prb} = rl; rhoUst{prb} = ru; frs_per_ang{prb} = fpa';

end

% selDir = selectivity index, 0 = cluster fires same for all directions and 1 = cluster only fired for one of the directions
Tclust.([EPOCH_NICKNAME, '_selDir']) = vertcat(sel_dir{:});

% prefDir = vector-averaged preferred direction in deg
Tclust.([EPOCH_NICKNAME, '_prefDir']) = vertcat(pref_dir{:});

% rhoLst = 5th CIs for each direction
Tclust.([EPOCH_NICKNAME, '_rhoLst']) = vertcat(rhoLst{:});

% rhoUst = 95th CIs for each direction
Tclust.([EPOCH_NICKNAME, '_rhoUst']) = vertcat(rhoUst{:});

% frs_perAng = firing rates for each direction
Tclust.([EPOCH_NICKNAME, '_frs_perAng']) = vertcat(frs_per_ang{:});

% mnFR_perAng = mean firing rate for each direction
Tclust.([EPOCH_NICKNAME, '_mnFR_perAng']) = cellfun(@(q) mean(q), Tclust.([EPOCH_NICKNAME, '_frs_perAng']), 'uni', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
