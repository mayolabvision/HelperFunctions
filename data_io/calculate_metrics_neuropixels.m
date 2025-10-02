function Snew = calculate_metrics_neuropixels(S)
    % VMI, if dirmem was ran
    fields = fieldnames(S);
    matchingFields1 = fields(contains(fields, {'dirmem', 'mdir'}, 'IgnoreCase', true));

    if ~isempty(matchingFields1)
        fprintf('\n~~CALCULATING MDIR METRICS~~\n');
        Tmdir = []; 
        for mm = 1:numel(matchingFields1)
            Tmdir = [Tmdir; S.(matchingFields1{mm}).tbl];
        end
        
        Tmdir = Tmdir(Tmdir.result=='CORRECT',:);


        for imec = 1:size(S.kilosort,1)
            % Percent trials fired on
            vals = cellfun(@(u) cellfun(@(q) numel(q), u, 'uni', 1)>0, Tmdir.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            S.kilosort(imec).clusters.ratio_mdir_trials = sum(vertcat(vals{:}),1)'./(height(Tmdir));
            
            % delay period FR
            delay_fr = cellfun(@(u,v,w) cellfun(@(q) (sum(q>=u(1) & q<v)./(v-u(1)))*1000, w, 'uni', 0), Tmdir.TARG_ON, num2cell(Tmdir.SACCADE), Tmdir.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            delay_fr = cell2mat(vertcat(delay_fr{:}));

            % dprime between evoked and spotantaneous activity
            vis_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u(1)+50) & q<(u(1)+250)), w, 'uni', 0), Tmdir.TARG_ON, Tmdir.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            vis_fr = cell2mat(vertcat(vis_fr{:}));

            sac_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u-75) & q<(u+75)), w, 'uni', 0), num2cell(Tmdir.SACCADE), Tmdir.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            sac_fr = cell2mat(vertcat(sac_fr{:}));

            spont_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u(1)-200) & q<(u(1)-50)), w, 'uni', 0), Tmdir.TARG_ON, Tmdir.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            spont_fr = cell2mat(vertcat(spont_fr{:}));

            dp_vis = (mean(vis_fr,1) - mean(spont_fr,1)) ./ (sqrt(0.5*(var(vis_fr,1) + var(spont_fr,1))));
            dp_sac = (mean(sac_fr,1) - mean(spont_fr,1)) ./ (sqrt(0.5*(var(sac_fr,1) + var(spont_fr,1))));

            S.kilosort(imec).clusters.dp_vis = dp_vis';
            S.kilosort(imec).clusters.dp_sac = dp_sac';

            %
            dirs = sort(unique(Tmdir.angle));
            fr_perDir = zeros(numel(dirs),size(delay_fr,2));
            for d = 1:numel(dirs)
                fr_perDir(d,:) = mean(delay_fr(Tmdir.angle==dirs(d),:),1);
            end

            S.kilosort(imec).clusters.mdir_delayFR_perDir = num2cell(fr_perDir');
            [ii,mm] = max(fr_perDir);
            S.kilosort(imec).clusters.mdir_delayFR_peakDirFR = ii';
            S.kilosort(imec).clusters.mdir_delayFR_peakDir = dirs(mm);

            % VISUAL
            [vis_sel_dir, vis_pref_dir, ~, ~, frs_perAng_vis] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[50,150],'ALIGN_TO','stim','IMEC',imec);
            S.kilosort(imec).clusters.vis_sel_dir = vis_sel_dir;
            S.kilosort(imec).clusters.vis_pref_dir = vis_pref_dir;

            % MOTOR
            [sac_sel_dir, sac_pref_dir, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','IMEC',imec);
            S.kilosort(imec).clusters.sac_sel_dir = sac_sel_dir;
            S.kilosort(imec).clusters.sac_pref_dir = sac_pref_dir;

            % VMI
            [VMI_per_unit,VMIdp_per_unit, visFR_per_unit, sacFR_per_unit] = deal(zeros(size(frs_perAng_sac,2),1));
            for unit = 1:size(frs_perAng_sac,2)
                visFR = frs_perAng_vis(:,unit);
                sacFR = frs_perAng_sac(:,unit);

                visFR_per_unit(unit) = mean(vertcat(visFR{:}));
                sacFR_per_unit(unit) = mean(vertcat(sacFR{:}));

                VMI_per_unit(unit) = (mean(vertcat(visFR{:})) - mean(vertcat(sacFR{:})))/(mean(vertcat(visFR{:})) + mean(vertcat(sacFR{:})));
                VMIdp_per_unit(unit) = (mean(vertcat(visFR{:}))-mean(vertcat(sacFR{:})))/sqrt(var(vertcat(visFR{:})) * var(vertcat(sacFR{:})));
            end
            S.kilosort(imec).clusters.vis_meanFR = visFR_per_unit;
            S.kilosort(imec).clusters.sac_meanFR = sacFR_per_unit;
            S.kilosort(imec).clusters.VMI = VMI_per_unit;
            S.kilosort(imec).clusters.VMIdp = VMIdp_per_unit;
        end
    end

    % PURSUIT
    matchingFields2 = fields(contains(fields, {'pursuit', 'purs'}, 'IgnoreCase', true));

    if ~isempty(matchingFields2)
        fprintf('\n~~CALCULATING PURS METRICS~~\n');
        Tpurs = []; 
        for mm = 1:numel(matchingFields2)
            Tpurs = [Tpurs; S.(matchingFields2{mm}).tbl];
        end
        Tpurs = Tpurs(Tpurs.result=='CORRECT' & Tpurs.jump==-1 & Tpurs.pursType=='pure' & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0),:);

        for imec = 1:size(S.kilosort,1)
            % Percent trials fired on
            vals = cellfun(@(u) cellfun(@(q) numel(q), u, 'uni', 1)>0, Tpurs.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            S.kilosort(imec).clusters.ratio_purs_trials = sum(vertcat(vals{:}),1)'./(height(Tpurs));

            % delay period FR
            targ_fr = cellfun(@(u,v,w) cellfun(@(q) (sum(q>=u & q<v)./(v-u(1)))*1000, w, 'uni', 0), num2cell(Tpurs.PURSUIT_TARG_ON), num2cell(Tpurs.PURSUIT_TARG_OFF), Tpurs.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            targ_fr = cell2mat(vertcat(targ_fr{:}));

            % dprime between evoked and spotantaneous activity
            pur_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u+50) & q<(u+250)), w, 'uni', 0), num2cell(Tpurs.PURSUIT_TARG_ON), Tpurs.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            pur_fr = cell2mat(vertcat(pur_fr{:}));

            spont_fr = cellfun(@(u,w) cellfun(@(q) sum(q>=(u-200) & q<(u-50)), w, 'uni', 0), num2cell(Tpurs.PURSUIT_TARG_ON), Tpurs.(sprintf('spiketimes_%d',imec)), 'uni', 0);
            spont_fr = cell2mat(vertcat(spont_fr{:}));

            dp_pur = (mean(pur_fr,1) - mean(spont_fr,1)) ./ (sqrt(0.5*(var(pur_fr,1) + var(spont_fr,1))));

            S.kilosort(imec).clusters.dp_pur = dp_pur';

            dirs = sort(unique(Tpurs.angle));
            fr_perDir = zeros(numel(dirs),size(targ_fr,2));
            for d = 1:numel(dirs)
                fr_perDir(d,:) = mean(targ_fr(Tpurs.angle==dirs(d),:),1);
            end

            S.kilosort(imec).clusters.purs_targFR_perDir = num2cell(fr_perDir');
            [ii,mm] = max(fr_perDir);
            S.kilosort(imec).clusters.purs_targFR_peakDirFR = ii';
            S.kilosort(imec).clusters.purs_targFR_peakDir = dirs(mm);

            % MOTOR (PURSUIT)
            [pur_sel_dir, pur_pref_dir, ~, ~, ~] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','IMEC',imec);
            S.kilosort(imec).clusters.pur_sel_dir = pur_sel_dir;
            S.kilosort(imec).clusters.pur_pref_dir = pur_pref_dir;
        end
    end

    % SPI, if dirmem and pursuit were ran
    if ~isempty(matchingFields1) & ~isempty(matchingFields2)
        for imec = 1:size(S.kilosort,1)
            % MOTOR (SACCADE)
            [~, ~, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','IMEC',imec);

            % MOTOR (PURSUIT)
            [~, ~, ~, ~, frs_perAng_pur] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','IMEC',imec);

            % VMI
            [SPI_per_unit,SPIdp_per_unit,purFR_per_unit] = deal(zeros(size(frs_perAng_sac,2),1));
            for unit = 1:size(frs_perAng_sac,2)
                sacFR = frs_perAng_sac(:,unit);
                purFR = frs_perAng_pur(:,unit);

                purFR_per_unit(unit) = mean(vertcat(purFR{:}));
            
                SPI_per_unit(unit) = (mean(vertcat(sacFR{:})) - mean(vertcat(purFR{:})))/(mean(vertcat(sacFR{:})) + mean(vertcat(purFR{:})));
                SPIdp_per_unit(unit) = (mean(vertcat(sacFR{:}))-mean(vertcat(purFR{:})))/sqrt(var(vertcat(sacFR{:})) * var(vertcat(purFR{:})));
            end
            S.kilosort(imec).clusters.pur_meanFR = purFR_per_unit;
            S.kilosort(imec).clusters.SPI = SPI_per_unit;
            S.kilosort(imec).clusters.SPIdp = SPIdp_per_unit;
        end
    end
    Snew = S;
end
