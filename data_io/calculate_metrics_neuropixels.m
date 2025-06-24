function Snew = calculate_metrics_neuropixels(S,trlAvg_frs_all)
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
            if ~isempty(trlAvg_frs_all{imec})
                % VISUAL
                [vis_sel_dir, vis_pref_dir, ~, ~, frs_perAng_vis] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[50,150],'ALIGN_TO','stim','IMEC',imec-1);
                S.kilosort(imec).clusters.vis_sel_dir = vis_sel_dir;
                S.kilosort(imec).clusters.vis_pref_dir = vis_pref_dir;

                % MOTOR
                [sac_sel_dir, sac_pref_dir, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','IMEC',imec-1);
                S.kilosort(imec).clusters.sac_sel_dir = sac_sel_dir;
                S.kilosort(imec).clusters.sac_pref_dir = sac_pref_dir;

                % VMI
                VMI_per_unit = zeros(size(frs_perAng_sac,2),1);
                for unit = 1:size(frs_perAng_sac,2)
                    visFR = frs_perAng_vis(:,unit);
                    sacFR = frs_perAng_sac(:,unit);
                
                    VMI_per_unit(unit) = (mean(vertcat(visFR{:})) - mean(vertcat(sacFR{:})))/(mean(vertcat(visFR{:})) + mean(vertcat(sacFR{:})));
                end
                S.kilosort(imec).clusters.VMI = VMI_per_unit;
            end
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
            if ~isempty(trlAvg_frs_all{imec})
                % MOTOR (PURSUIT)
                [pur_sel_dir, pur_pref_dir, ~, ~, ~] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','IMEC',imec-1);
                S.kilosort(imec).clusters.pur_sel_dir = pur_sel_dir;
                S.kilosort(imec).clusters.pur_pref_dir = pur_pref_dir;
            end
        end
    end

    % SPII, if dirmem and pursuit were ran
    if ~isempty(matchingFields1) & ~isempty(matchingFields2)
        for imec = 1:size(S.kilosort,1)
            if ~isempty(trlAvg_frs_all{imec})
                % MOTOR (SACCADE)
                [~, ~, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','IMEC',imec-1);

                % MOTOR (PURSUIT)
                [~, ~, ~, ~, frs_perAng_pur] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','IMEC',imec-1);

                % VMI
                SPI_per_unit = zeros(size(frs_perAng_sac,2),1);
                for unit = 1:size(frs_perAng_sac,2)
                    sacFR = frs_perAng_sac(:,unit);
                    purFR = frs_perAng_pur(:,unit);
                
                    SPI_per_unit(unit) = (mean(vertcat(sacFR{:})) - mean(vertcat(purFR{:})))/(mean(vertcat(sacFR{:})) + mean(vertcat(purFR{:})));
                end
                S.kilosort(imec).clusters.SPI = SPI_per_unit;
            end
        end
    end
    Snew = S;
end