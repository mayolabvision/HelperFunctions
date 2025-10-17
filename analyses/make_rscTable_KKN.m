function make_rscTable_KKN(data_path)

load(data_path, 'C', 'Tmdir', 'Tpurs', 'Trfmp');

fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

Tmdir.sess_name = categorical(regexprep(string(Tmdir.sess_name),'_g0$',''));
Tpurs.sess_name = categorical(regexprep(string(Tpurs.sess_name),'_g0$',''));
Trfmp.sess_name = categorical(regexprep(string(Trfmp.sess_name),'_g0$',''));
C.sess_name = categorical(regexprep(string(C.sess_name),'_g0$',''));

% Removing any units with NaN or Inf in selectivity metrics
C1 = C(~isnan(C.VMI) & ~isinf(C.VMI) & ~isnan(C.VMIdp) & ~isinf(C.VMIdp) & ~isnan(C.SPI) & ...
       ~isinf(C.SPI) & ~isnan(C.SPIdp) & ~isinf(C.SPIdp) & ~isnan(C.SPI) & ~isinf(C.SPI) & ...
       C.vis_sel_dir>0 & C.vis_sel_dir<1 & C.sac_sel_dir>0 & C.sac_sel_dir<1 & C.pur_sel_dir>0 & C.pur_sel_dir<1,:);

% Removing any units that don't fire at least once on most of the trials
PR = 0.95;
C2 = C1(C1.ratio_mdir_trials >= PR & C1.ratio_purs_trials >=PR,:);

% Minimum dp between evoked and spontaneous activity
DP = 0.3;
C3 = C2(C2.dp_vis >= DP | C2.dp_sac >= DP | C2.dp_pur >= DP,:);

% Minimum FR
FR = 1;
C4 = C3(C3.mdir_delayFR_peakDirFR >= FR | C3.purs_targFR_peakDirFR >= FR,:);

% Minimum SNR
SNR = 2;
C5 = C4(C4.snr>=SNR,:);

CC = C5;

nBins = 5;
CC.VMI_quantile = discretize(CC.VMI,  quantile(CC.VMI, linspace(0, 1, nBins+1)));
CC.VMI_quantile(CC.VMI == max(CC.VMI)) = nBins;

CC.SPI_quantile = discretize(CC.SPI, quantile(CC.SPI, linspace(0, 1, nBins+1)));
CC.SPI_quantile(CC.SPI == max(CC.SPI)) = nBins;

sess = unique(CC.sess_name);
pairs_all = cell(length(sess),1);

for s = 1:length(sess)
    fprintf('%d of %d\n', s, numel(sess));
    this_mdir = Tmdir(Tmdir.sess_name==sess(s),:);
    this_purs = Tpurs(Tpurs.sess_name==sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);

    vis_frs_l = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w(1)+50 & u<w(1)+250))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_1, this_mdir.TARG_ON, 'uni', 0);
    vis_frs_r = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w(1)+50 & u<w(1)+250))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_2, this_mdir.TARG_ON, 'uni', 0);

    sac_frs_l = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_1, num2cell(this_mdir.SACCADE), 'uni', 0);
    sac_frs_r = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_mdir.spiketimes_2, num2cell(this_mdir.SACCADE), 'uni', 0);

    pur_frs_l = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_purs.spiketimes_1, num2cell(this_purs.pursuitOnset), 'uni', 0);
    pur_frs_r = cellfun(@(q,w) cellfun(@(u) ((sum(u>=w-100 & u<w+100))/(200))*1000, q, 'uni', 0), this_purs.spiketimes_2, num2cell(this_purs.pursuitOnset), 'uni', 0);
    
    pairs_probes = cell(4,1);
    for p1 = 1:2
        if p1==1
            vis_frs_p1 = vis_frs_l;
            sac_frs_p1 = sac_frs_l;
            pur_frs_p1 = pur_frs_l;
        else
            vis_frs_p1 = vis_frs_r;
            sac_frs_p1 = sac_frs_r;
            pur_frs_p1 = pur_frs_r;
        end

        for p2 = 1:2
            if p2==1
                vis_frs_p2 = vis_frs_l;
                sac_frs_p2 = sac_frs_l;
                pur_frs_p2 = pur_frs_l;
            else
                vis_frs_p2 = vis_frs_r;
                sac_frs_p2 = sac_frs_r;
                pur_frs_p2 = pur_frs_r;
            end

            fprintf('%d - %d\n', p1, p2);

            p1_units = CC(CC.sess_name==sess(s) & CC.probe_index==p1,:);
            p2_units = CC(CC.sess_name==sess(s) & CC.probe_index==p2,:);
        
            if p1==p2
                pairs = cell((height(p1_units)*height(p2_units)) - height(p1_units),29);
            else
                pairs = cell(height(p1_units)*height(p2_units),29);
            end

            rr = 1;
            for n1 = 1:height(p1_units)
                for n2 = 1:height(p2_units)  
                    if ~(p1==p2 && n1==n2)
                        n1_clust = p1_units.cluster_id(n1);
                        n2_clust = p2_units.cluster_id(n2);
    
                        % SACCADE TASK
                        TT = this_mdir; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        angs = sort(unique(TT.angle));
                        amps = sort(unique(TT.distance)); 
    
                        %------------------------------ VISUAL ------------------------------%
                        frs_p1 = vis_frs_p1; frs_p2 = vis_frs_p2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if p1_units.vis_sel_dir(n1) > p2_units.vis_sel_dir(n2) %%%%%%%%%%%%%%%%%%
                            bestDir = p1_units.vis_pref_dir(n1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            bestDir = p2_units.vis_pref_dir(n2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
    
                        relAngs = mod(angs - bestDir + 180, 360) - 180;
                        [~, sortIdx] = sort(relAngs);
    
                        [rsc_perDir, fr_perDir] = deal(zeros(1, numel(angs)));
                        n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                        for d = 1:numel(angs)
                            thisAng = angs(d);
    
                            n12_fr = []; n1_zfr = []; n2_zfr = [];
                            for a = 1:numel(amps)
                                n1_fr = cellfun(@(q) q{n1_clust+1}, frs_p1(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
                                n2_fr = cellfun(@(q) q{n2_clust+1}, frs_p2(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
    
                                n12_fr = [n12_fr; n1_fr; n2_fr];
    
                                n1_zfr = [n1_zfr; zscore(n1_fr)];
                                n2_zfr = [n2_zfr; zscore(n2_fr)];
                            end
                            n1_zfr_all = [n1_zfr_all; n1_zfr];
                            n2_zfr_all = [n2_zfr_all; n2_zfr];
                            n12_fr_all = [n12_fr_all; n12_fr];
    
                            [rho,~] = corr(n1_zfr, n2_zfr);
    
                            rsc_perDir(sortIdx(d)) = rho;
                            fr_perDir(sortIdx(d)) = mean(n12_fr);
                        end
    
                        [rsig,~] = corr(cellfun(@(q) q{n1_clust+1}, frs_p1, 'uni', 1), cellfun(@(q) q{n2_clust+1}, frs_p2, 'uni', 1));
                        [rsc,~] = corr(n1_zfr_all, n2_zfr_all);
                        mnFR = mean(n12_fr_all);
    
                        rsig_vis = rsig; rsc_vis = rsc; mnFR_vis = mnFR; rsc_perDir_vis = rsc_perDir; fr_perDir_vis = fr_perDir; 
    
                        %------------------------------ SACCADE ------------------------------%
                        frs_p1 = sac_frs_p1; frs_p2 = sac_frs_p2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if p1_units.sac_sel_dir(n1) > p2_units.sac_sel_dir(n2) %%%%%%%%%%%%%%%%%%
                            bestDir = p1_units.sac_pref_dir(n1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else
                            bestDir = p2_units.sac_pref_dir(n2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
    
                        relAngs = mod(angs - bestDir + 180, 360) - 180;
                        [~, sortIdx] = sort(relAngs);
    
                        [rsc_perDir, fr_perDir] = deal(zeros(1, numel(angs)));
                        n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                        for d = 1:numel(angs)
                            thisAng = angs(d);
    
                            n12_fr = []; n1_zfr = []; n2_zfr = [];
                            for a = 1:numel(amps)
                                n1_fr = cellfun(@(q) q{n1_clust+1}, frs_p1(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
                                n2_fr = cellfun(@(q) q{n2_clust+1}, frs_p2(TT.angle==angs(d) & TT.distance==amps(a)), 'uni', 1);
    
                                n12_fr = [n12_fr; n1_fr; n2_fr];
    
                                n1_zfr = [n1_zfr; zscore(n1_fr)];
                                n2_zfr = [n2_zfr; zscore(n2_fr)];
                            end
                            n1_zfr_all = [n1_zfr_all; n1_zfr];
                            n2_zfr_all = [n2_zfr_all; n2_zfr];
                            n12_fr_all = [n12_fr_all; n12_fr];
    
                            [rho,~] = corr(n1_zfr, n2_zfr);
    
                            rsc_perDir(sortIdx(d)) = rho;
                            fr_perDir(sortIdx(d)) = mean(n12_fr);
                        end
    
                        [rsig,~] = corr(cellfun(@(q) q{n1_clust+1}, frs_p1, 'uni', 1), cellfun(@(q) q{n2_clust+1}, frs_p2, 'uni', 1));
                        [rsc,~] = corr(n1_zfr_all, n2_zfr_all);
                        mnFR = mean(n12_fr_all);
    
                        rsig_sac = rsig; rsc_sac = rsc; mnFR_sac = mnFR; rsc_perDir_sac = rsc_perDir; fr_perDir_sac = fr_perDir;
    
                        %------------------------------ PURSUIT ------------------------------%
                        TT = this_purs;
                        angs = sort(unique(TT.angle));
                        spes = sort(unique(TT.pursuitSpeed));
    
                        frs_p1 = pur_frs_p1; frs_p2 = pur_frs_p2; 
                        if p1_units.pur_sel_dir(n1) > p2_units.pur_sel_dir(n2) 
                            bestDir = p1_units.pur_pref_dir(n1);
                        else
                            bestDir = p2_units.pur_pref_dir(n2); 
                        end
    
                        relAngs = mod(angs - bestDir + 180, 360) - 180;
                        [~, sortIdx] = sort(relAngs);
    
                        [rsc_perDir, fr_perDir] = deal(zeros(1, numel(angs)));
                        n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                        for d = 1:numel(angs)
                            thisAng = angs(d);
    
                            n12_fr = []; n1_zfr = []; n2_zfr = [];
                            for a = 1:numel(spes)
                                n1_fr = cellfun(@(q) q{n1_clust+1}, frs_p1(TT.angle==angs(d) & TT.pursuitSpeed==spes(a)), 'uni', 1);
                                n2_fr = cellfun(@(q) q{n2_clust+1}, frs_p2(TT.angle==angs(d) & TT.pursuitSpeed==spes(a)), 'uni', 1);
    
                                n12_fr = [n12_fr; n1_fr; n2_fr];
    
                                n1_zfr = [n1_zfr; zscore(n1_fr)];
                                n2_zfr = [n2_zfr; zscore(n2_fr)];
                            end
                            n1_zfr_all = [n1_zfr_all; n1_zfr];
                            n2_zfr_all = [n2_zfr_all; n2_zfr];
                            n12_fr_all = [n12_fr_all; n12_fr];
    
                            [rho,~] = corr(n1_zfr, n2_zfr);
    
                            rsc_perDir(sortIdx(d)) = rho;
                            fr_perDir(sortIdx(d)) = mean(n12_fr);
                        end
    
                        [rsig,~] = corr(cellfun(@(q) q{n1_clust+1}, frs_p1, 'uni', 1), cellfun(@(q) q{n2_clust+1}, frs_p2, 'uni', 1));
                        [rsc,~] = corr(n1_zfr_all, n2_zfr_all);
                        mnFR = mean(n12_fr_all);
    
                        rsig_pur = rsig; rsc_pur = rsc; mnFR_pur = mnFR; rsc_perDir_pur = rsc_perDir; fr_perDir_pur = fr_perDir;
    
                        %---------------------------------------------------------------------------------------------------------
                        
                        pairs(rr,:) = {p1_units.monkey(n1), p1_units.sess_name(n1) ...
                                p1, p2, n1_clust, n2_clust, ...
                                p1_units.VMI_quantile(n1), p2_units.VMI_quantile(n2), ...
                                p1_units.SPI_quantile(n1), p2_units.SPI_quantile(n2), ...
                                norm(p1_units.unit_locations{n1}(1:2) - p2_units.unit_locations{n2}(1:2)), ...
                                abs(mod(p1_units.vis_pref_dir(n1) - p2_units.vis_pref_dir(n2) + 180, 360) - 180), ...
                                abs(mod(p1_units.sac_pref_dir(n1) - p2_units.sac_pref_dir(n2) + 180, 360) - 180), ...
                                abs(mod(p1_units.pur_pref_dir(n1) - p2_units.pur_pref_dir(n2) + 180, 360) - 180), ...
                                mnFR_vis, rsc_vis, rsig_vis, rsc_perDir_vis, fr_perDir_vis, ...
                                mnFR_sac, rsc_sac, rsig_sac, rsc_perDir_sac, fr_perDir_sac, ...
                                mnFR_pur, rsc_pur, rsig_pur, rsc_perDir_pur, fr_perDir_pur
                                };
                        rr = rr + 1;
                    end
                end
            end   
            pairs_probes{(p1-1)*2 + p2} = pairs;
        end
    end
    pairs_all{s} = pairs_probes;
end

pairs_all2 = vertcat(pairs_all{:});
pairs_all3 = vertcat(pairs_all2{:});
pairs_tbl = cell2table(pairs_all3,'VariableNames',{'monkey','sess_name','n1_probe','n2_probe','n1_cluster','n2_cluster','n1_VMI','n2_VMI','n1_SPI','n2_SPI','depth_diff','vis_prefDir_diff','sac_prefDir_diff','pur_prefDir_diff','mnFR_vis','rsc_vis','rsig_vis','rsc_perDir_vis','fr_perDir_vis','mnFR_sac','rsc_sac','rsig_sac','rsc_perDir_sac','fr_perDir_sac','mnFR_pur','rsc_pur','rsig_pur','rsc_perDir_pur','fr_perDir_pur'});

save('/ix1/pmayo/lab_NHPdata/pairs_table.mat', 'pairs_tbl');

end
