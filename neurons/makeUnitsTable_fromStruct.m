function [unitsTbl] = makeUnitsTable_fromStruct(exp,trialTbl,postint)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    unitnames = exp.info.channels;
    snrs = exp.info.SNRs;

    % If >24 then MT, <24 then FEF
    unitNum     =  cellfun(@(y) sscanf(y,'unit%d')>24, unitnames,'uni',1);
    brainareas  =  cell(length(unitNum),1);
    brainareas(unitNum==1,1)  =  {'MT'};
    brainareas(unitNum==0,1)  =  {'FEF'};

    dirsdeg = sort(unique(trialTbl.Direction));
    speeds = sort(unique(trialTbl.Speed));

    spks = cellfun(@(q) struct2cell(q), {exp.dataMaestroPlx.units}.','uni', 0);
    spks = horzcat(spks{:});

    % pull out pure pursuit trials only
    unitsTbl = cell(length(unitnames),15);
    for u=1:length(unitnames)
        spkcnts = cell(1,length(dirsdeg));
        for d=1:length(dirsdeg)
            spkcnts{d} = cellfun(@(q,s) (sum(q>=s & q<s+postint)), spks(u,trialTbl.Direction==dirsdeg(d)), num2cell(trialTbl.TargetMotionOnset(trialTbl.Direction==dirsdeg(d))'),'uni', 1); 
        end
        spkcnts_bySpeeds = cell(1,length(speeds));
        for s=1:length(speeds)
            spkcnts_bySpeeds{s} = cellfun(@(q,s) (sum(q>=s & q<s+postint)), spks(u,trialTbl.Speed==speeds(s)), num2cell(trialTbl.TargetMotionOnset(trialTbl.Speed==speeds(s))'),'uni', 1); 
        end

        mnFRByDir= cellfun(@nanmean, spkcnts);
        varFRByDir= cellfun(@(q) var(q,'omitnan'), spkcnts, 'uni', 1);
        [mnFR_bestDir,m] = max(mnFRByDir);
        bestDir = dirsdeg(m);
        varFR_bestDir = varFRByDir(m);

        mnFRBySp= cellfun(@nanmean, spkcnts_bySpeeds);

        if length(dirsdeg)==4
            if m==1
                n = 3;
            elseif m==2
                n = 4;
            elseif m==3
                n = 1;
            elseif m==4
                n = 2;
            end
        elseif length(dirsdeg)==2
            if m==1
                n = 2;
            elseif m==2
                n = 1;
            end
        end

        nullDir = dirsdeg(n);
        
        depthMod = (max(mnFRByDir) - min(mnFRByDir))./(max(mnFRByDir) + min(mnFRByDir)); % contrast ratio
        DI = 1 - (mnFRByDir(n)/mnFRByDir(m));

        [seldir, prefdirfit] = tuningbias(dirsdeg, mnFRByDir); % Matt Smith's 'tuningbias' function to estimate best direction; prefdir must be in DEGREES
   
        % Selectivity Index (SI):
        if length(dirsdeg)==4
            resp_ang = num2cell([horzcat(spkcnts{:})' [repmat(dirsdeg(1),1,length(spkcnts{1})),repmat(dirsdeg(2),1,length(spkcnts{2})),repmat(dirsdeg(3),1,length(spkcnts{3})),repmat(dirsdeg(4),1,length(spkcnts{4}))]'],2);
        elseif length(dirsdeg)==2
            resp_ang = num2cell([horzcat(spkcnts{:})' [repmat(dirsdeg(1),1,length(spkcnts{1})),repmat(dirsdeg(2),1,length(spkcnts{2}))]'],2);
        end
        SI = sqrt(((sum(cellfun(@(y) y(1)*sin(deg2rad(y(2))), resp_ang, 'uni', 1)))^2) + ((sum(cellfun(@(z) z(1)*cos(deg2rad(z(2))), resp_ang, 'uni', 1)))^2))/sum(cellfun(@(q) q(1), resp_ang,'uni',1));
       
        % Generate randomized index of stimrate values, WITH REPLACEMENT
        shuffles = 1000;
        permFR = nan(shuffles,size(spkcnts, 1) ); % pre-allocate
        for c = 1:size(spkcnts, 1) % for each condition
            for sh=1:shuffles % for each shuffle
                permFR(sh,c) = mean(randsample(cell2mat(spkcnts), size(spkcnts{c}, 2), 1 )); % Select X number of spks/trls WITH replacement
            end
        end
       
        permFR = sort(permFR);
        rhoLst = quantile(permFR, 0.025); % 95% Lower confidence interval
        rhoUst = quantile(permFR, 0.975); % 95% Upper confidence interval
        
        % if any of the FR are greater rhoUst or less than Lst then set signif flag
        if any (mnFRByDir > rhoUst) ||  any (mnFRByDir < rhoLst)
            signiffl = 1;
        else
            signiffl = 0;
        end
        
        %%%%%%%%%%%%%
        unitsTbl(u,:) = {unitnames{u}, [exp.info.expName,'_',unitnames{u}], brainareas{u}, snrs(u), bestDir, nullDir, prefdirfit, mnFR_bestDir, varFR_bestDir, depthMod, seldir, DI, SI, signiffl, {spks(u,:)}};
    end
end