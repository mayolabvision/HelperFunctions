function [sel_dir,pref_dir,rhoLst,rhoUst,frs_per_ang] = calculate_direction_tuning_from_tbl(tbl,varargin)
    p = inputParser;
    addRequired(p, 'tbl', @istable);
    addParameter(p, 'FR_WIN', [50,150], @isnumeric);
    addParameter(p, 'ALIGN_TO', 'stim', @ischar);
    addParameter(p, 'IMEC', 0, @isnumeric);
    
    parse(p, tbl, varargin{:});
    FR_WIN = p.Results.FR_WIN;
    ALIGN_TO = p.Results.ALIGN_TO;
    IMEC = p.Results.IMEC;

    theta = sort(unique(tbl.angle))';

    spikes = tbl.(sprintf('spiketimes_%d',IMEC));
    frs_per_ang = cell(length(theta),length(spikes{1}));
    for a = 1:length(theta)
        this_ang = tbl(tbl.angle==theta(a),:);
        if isequal(ALIGN_TO,'stim')
            FR = cellfun(@(w,v) cellfun(@(q) sum(q>=(v(1)+FR_WIN(1)) & q<(v(1)+FR_WIN(2)))/((FR_WIN(2)-FR_WIN(1))/1000), w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), this_ang.TARG_ON, 'uni', 0);
        elseif isequal(ALIGN_TO,'sacc')
            FR = cellfun(@(w,v) cellfun(@(q) sum(q>=(v+FR_WIN(1)) & q<(v+FR_WIN(2)))/((FR_WIN(2)-FR_WIN(1))/1000), w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), num2cell(this_ang.SACCADE), 'uni', 0);
        elseif isequal(ALIGN_TO,'targ')
            FR = cellfun(@(w,v) cellfun(@(q) sum(q>=(v+FR_WIN(1)) & q<(v+FR_WIN(2)))/((FR_WIN(2)-FR_WIN(1))/1000), w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), num2cell(this_ang.PURSUIT_TARG_ON), 'uni', 0);
        elseif isequal(ALIGN_TO,'purs')
            FR = cellfun(@(w,v) cellfun(@(q) sum(q>=(v+FR_WIN(1)) & q<(v+FR_WIN(2)))/((FR_WIN(2)-FR_WIN(1))/1000), w, 'uni', 0), this_ang.(sprintf('spiketimes_%d',IMEC)), num2cell(this_ang.pursuitOnset), 'uni', 0);
        end
        frs_per_ang(a,:) = num2cell(cell2mat(vertcat(FR{:})),1);
    end
    
    maxLength = max(cellfun(@numel, frs_per_ang)); maxLength = max(maxLength);
    [sel_dir, pref_dir] = deal(zeros(size(frs_per_ang,2),1));
    for unit = 1:size(frs_per_ang,2)
        frs_perAng2 = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_per_ang(:,unit), 'UniformOutput', false);
                        
        stimrate = vertcat(frs_perAng2{:})';
        
        % Generate randomized index of stimrate values, WITH REPLACEMENT
        shuffles = 1000;
        rhoPst = [];
        
        for sh=1:shuffles
            randind=randi( (size(stimrate,1)*size(stimrate,2)), size(stimrate,1), size(stimrate,2) );
            permutedStimrate = stimrate(randind);
            rhoPst = [rhoPst; mean(permutedStimrate, 'omitnan')];
        end
        
        sorted_rhoPst=sort(rhoPst);
        rhoLst = sorted_rhoPst(shuffles*.05,:); % 95% lower confidence interval
        rhoUst = sorted_rhoPst(shuffles-(shuffles*.05),:); % 95% upper confidence interval
        
        % calculate tuning preferences
        [visds, visdp] = tuningbias(theta,mean(stimrate,'omitnan'));
        sel_dir(unit) = visds; 
        pref_dir(unit) = visdp;
    end
end