function [sel_dir,pref_dir,rhoLst,rhoUst,frs_per_ang,pval_dir] = calculate_tuning(tbl, varargin)
    p = inputParser;
    addRequired(p, 'tbl', @istable);
    addParameter(p, 'PROBE_INDEX', [], @isnumeric);
    addParameter(p, 'FR_WIN', [50,150], @isnumeric);
    addParameter(p, 'ALIGN_TO', 'TARG_ON', @ischar);
    addParameter(p, 'OPT_LAT', [], @isnumeric);
    addParameter(p, 'nShuffles', 1000, @isnumeric);
    addParameter(p, 'WITH_BONF', false, @islogical);
    addParameter(p, 'WITH_MAXSTAT', true, @islogical);
    
    parse(p, tbl, varargin{:});
    PROBE_INDEX = p.Results.PROBE_INDEX;
    FR_WIN = p.Results.FR_WIN;
    ALIGN_TO = p.Results.ALIGN_TO;
    OPT_LAT = p.Results.OPT_LAT;
    nShuffles = p.Results.nShuffles;
    WITH_BONF = p.Results.WITH_BONF;
    WITH_MAXSTAT = p.Results.WITH_MAXSTAT;

    if isnumeric(tbl.(ALIGN_TO))
        tbl.(ALIGN_TO) = num2cell(tbl.(ALIGN_TO));
    end

    if ~isempty(PROBE_INDEX)
        prb_name = sprintf('spiketimes_%d',PROBE_INDEX);
    else
        prb_name = 'spiketimes';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract firing rates per direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = sort(unique(tbl.angle))'; % all possible target angles
    frs_per_ang = cell(length(theta),length(tbl.(prb_name){1}));

    for a = 1:length(theta)
        this_ang = tbl(tbl.angle==theta(a),:);

        if isempty(OPT_LAT)
            FR = cellfun(@(w,v) cellfun(@(q) ...
                sum(q >= (v(1)+FR_WIN(1)) & q < (v(1)+FR_WIN(2))) / ...
                ((FR_WIN(2)-FR_WIN(1))/1000), ...
                w, 'uni', 0), this_ang.(prb_name), ...
                this_ang.(ALIGN_TO), 'uni', 0);
        else
            FR = cellfun(@(w,v) cellfun(@(q,z) ...
                sum(q >= (v(1)+z+FR_WIN(1)) & q < (v(1)+z+FR_WIN(2))) / ...
                ((FR_WIN(2)-FR_WIN(1))/1000), ...
                w, num2cell(OPT_LAT'), 'uni', 0), this_ang.(prb_name), ...
                this_ang.(ALIGN_TO), 'uni', 0);
        end
           
        frs_per_ang(a,:) = num2cell(cell2mat(vertcat(FR{:})),1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute tuning and perform shuffle-based significance test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maxLength = max(cellfun(@numel, frs_per_ang)); maxLength = max(maxLength);

    [sel_dir, pref_dir, pval_dir] = deal(zeros(size(frs_per_ang,2),1));
    [rhoLst,rhoUst] = deal(cell(size(frs_per_ang,2),1));
    
    for unit = 1:size(frs_per_ang,2)
        frs_perAng2 = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_per_ang(:,unit), 'uni', 0);
                        
        %--------------- For a single unit -------------%
        % stimrate: trials × directions
        stimrate = vertcat(frs_perAng2{:})'; 
       
        % --- True direction tuning ---
        [ds, dp] = tuningbias(theta,mean(stimrate,'omitnan'));
        sel_dir(unit) = ds; 
        pref_dir(unit) = dp;
    
        % Observed mean firing rate per direction
        mnFR_obs = mean(stimrate,'omitnan');  
    
        % Shuffle null
        mnFR_null = zeros(nShuffles,numel(theta));
    
        % If using max-stat, also track scalar null
        if WITH_MAXSTAT
            stat_null = zeros(nShuffles,1);
        end
    
        for sh = 1:nShuffles
            shuffled = stimrate(randperm(numel(stimrate)));
            shuffled = reshape(shuffled, size(stimrate));
    
            mn = mean(shuffled,'omitnan');
            mnFR_null(sh,:) = mn;
    
            if WITH_MAXSTAT
                % --- Choose ONE of these stats ---
                % Option A: range (recommended)
                stat_null(sh) = max(mn) - min(mn);
    
                % Option B: max abs (alternative)
                % stat_null(sh) = max(abs(mn));
            end
        end
    
        % ================================
        % --- SIGNIFICANCE TEST OPTIONS ---
        % ================================
        if WITH_MAXSTAT
            % --- Observed stat ---
            stat_obs = max(mnFR_obs) - min(mnFR_obs);
            % stat_obs = max(abs(mnFR_obs)); % alternative
    
            % --- p-value ---
            pval_dir(unit) = double(mean(stat_null >= stat_obs) < 0.05);
    
            % Not used in this mode
            rhoLst{unit} = [];
            rhoUst{unit} = [];
    
        else
            if WITH_BONF
                alpha = 0.05;
                nDir = numel(theta);
                alpha_corr = alpha / nDir;
    
                lower_pct = 100 * (alpha_corr / 2);
                upper_pct = 100 * (1 - alpha_corr / 2);
    
                rhoLst{unit} = prctile(mnFR_null, lower_pct);
                rhoUst{unit} = prctile(mnFR_null, upper_pct);
    
            else
                % Per-direction CI (uncorrected)
                rhoLst{unit} = prctile(mnFR_null, 2.5);
                rhoUst{unit} = prctile(mnFR_null, 97.5);
            end
    
            % Any-direction logic
            pval_dir(unit) = double(any(mnFR_obs < rhoLst{unit} | mnFR_obs > rhoUst{unit}));
        end
    end
    
end
