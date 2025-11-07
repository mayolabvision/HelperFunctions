function [sel_dir,pref_dir,rhoLst,rhoUst,frs_per_ang, pval_dir, sig_dir] = calculate_direction_tuning_from_tbl(tbl,varargin)
    p = inputParser;
    addRequired(p, 'tbl', @istable);
    addParameter(p, 'FR_WIN', [50,150], @isnumeric);
    addParameter(p, 'ALIGN_TO', 'stim', @ischar);
    addParameter(p, 'PROBE', 0, @isnumeric);
    addParameter(p, 'nShuffles', 1000, @isnumeric);
    addParameter(p, 'alpha', 0.05, @isnumeric);
    
    parse(p, tbl, varargin{:});
    FR_WIN = p.Results.FR_WIN;
    ALIGN_TO = p.Results.ALIGN_TO;
    PROBE = p.Results.PROBE;
    nShuffles = p.Results.nShuffles;
    alpha = p.Results.alpha;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract firing rates per direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = sort(unique(tbl.angle))'; % all possible target angles
    spikes = tbl.(sprintf('spiketimes_%d',PROBE));
    frs_per_ang = cell(length(theta),length(spikes{1}));

    for a = 1:length(theta)
        this_ang = tbl(tbl.angle==theta(a),:);
        switch ALIGN_TO
            case 'stim'
                FR = cellfun(@(w,v) cellfun(@(q) ...
                    sum(q >= (v(1)+FR_WIN(1)) & q < (v(1)+FR_WIN(2))) / ...
                    ((FR_WIN(2)-FR_WIN(1))/1000), ...
                    w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), ...
                    this_ang.TARG_ON, 'uni', 0);
            case 'sacc'
                FR = cellfun(@(w,v) cellfun(@(q) ...
                    sum(q >= (v+FR_WIN(1)) & q < (v+FR_WIN(2))) / ...
                    ((FR_WIN(2)-FR_WIN(1))/1000), ...
                    w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), ...
                    num2cell(this_ang.saccadeOnset), 'uni', 0);
            case 'targ'
                FR = cellfun(@(w,v) cellfun(@(q) ...
                    sum(q >= (v+FR_WIN(1)) & q < (v+FR_WIN(2))) / ...
                    ((FR_WIN(2)-FR_WIN(1))/1000), ...
                    w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), ...
                    num2cell(this_ang.PURSUIT_TARG_ON), 'uni', 0);
            case 'purs'
                FR = cellfun(@(w,v) cellfun(@(q) ...
                    sum(q >= (v+FR_WIN(1)) & q < (v+FR_WIN(2))) / ...
                    ((FR_WIN(2)-FR_WIN(1))/1000), ...
                    w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), ...
                    num2cell(this_ang.pursuitOnset), 'uni', 0);
        end
        frs_per_ang(a,:) = num2cell(cell2mat(vertcat(FR{:})),1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute tuning and perform shuffle-based significance test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxLength = max(cellfun(@numel, frs_per_ang)); maxLength = max(maxLength);

    [sel_dir, pref_dir, rhoLst, rhoUst, pval_dir] = deal(zeros(size(frs_per_ang,2),1));
    for unit = 1:size(frs_per_ang,2)
        frs_perAng2 = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_per_ang(:,unit), 'uni', 0);
                        
        stimrate = vertcat(frs_perAng2{:})'; % trials × directions

        % --- True direction tuning ---
        [visds, visdp] = tuningbias(theta,mean(stimrate,'omitnan'));
        sel_dir(unit) = visds; 
        pref_dir(unit) = visdp;

        % --- Shuffle test ---
        rho_shuff = zeros(nShuffles,1);
        for sh = 1:nShuffles
            shuffled = stimrate(randperm(numel(stimrate)));
            shuffled = reshape(shuffled, size(stimrate));
            rho_shuff(sh) = tuningbias(theta, mean(shuffled, 'omitnan'));
        end

        % --- Compute CI bounds ---
        sorted_rho = sort(rho_shuff);
        rhoLst(unit) = sorted_rho(round(nShuffles*alpha/2));
        rhoUst(unit) = sorted_rho(round(nShuffles*(1-alpha/2)));

        % --- Compute p-value ---
        pval_dir(unit) = mean(rho_shuff >= true_rho);
    end

    sig_dir = pval_dir < alpha;
    
end