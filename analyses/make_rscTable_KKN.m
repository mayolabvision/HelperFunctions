function make_rscTable_KKN(data_path)

addpath(genpath('/ihome/pmayo/knoneman/HelperFunctions'))
load(data_path, 'C', 'Tmdir', 'Tpurs', 'Trfmp');

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

sess = unique(CC.sess_name);
pairs_all = cell(length(sess),1);

for s = 1:length(sess)
    fprintf('%d of %d\n', s, numel(sess));
    this_mdir = Tmdir(Tmdir.sess_name==sess(s),:);
    this_purs = Tpurs(Tpurs.sess_name==sess(s) & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0) & Tpurs.pursType=='pure',:);

    pairs_probes = cell(4,1);
    for p1 = 1:2
        tStart = tic;  % start timer
        for p2 = 1:2
            fprintf('%d - %d\n', p1, p2);

            p1_units = CC(CC.sess_name==sess(s) & CC.probe_index==p1,:);
            p2_units = CC(CC.sess_name==sess(s) & CC.probe_index==p2,:);
        
            pairs = cell(height(p1_units)*height(p2_units),24);
            rr = 1;
            for n1 = 1:height(p1_units)
                for n2 = 1:height(p2_units)
                    % MDIR
                    angs = sort(unique(this_mdir.angle));
                    amps = sort(unique(this_mdir.distance));

                    if 0.5*(p1_units.vis_sel_dir(n1)+p1_units.sac_sel_dir(n1)) > 0.5*(p2_units.vis_sel_dir(n2)+p2_units.sac_sel_dir(n2))
                        bestDir = p1_units.mdir_delayFR_peakDir(n1);
                    else
                        bestDir = p2_units.mdir_delayFR_peakDir(n2);
                    end

                    relAngs = mod(angs - bestDir + 180, 360) - 180;
                    [~, sortIdx] = sort(relAngs);

                    [rsc_perDir_mdir, fr_perDir_mdir] = deal(zeros(1, numel(angs)));
                    n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                    for d = 1:numel(angs)
                        thisAng = angs(d);

                        n12_fr = []; n1_zfr = []; n2_zfr = [];
                        for a = 1:numel(amps)
                            n12_fr = [n1_fr; [cellfun(@(q,w,v) ((sum(q{n1}>w(1) & q{n1}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p1))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1);cellfun(@(q,w,v) ((sum(q{n2}>w(1) & q{n2}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p2))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1)]];
                            
                            n1_zfr = [n1_zfr; zscore(cellfun(@(q,w,v) ((sum(q{n1}>w(1) & q{n1}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p1))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1))];
                            n2_zfr = [n2_zfr; zscore(cellfun(@(q,w,v) ((sum(q{n2}>w(1) & q{n2}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p2))(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), this_mdir.TARG_ON(this_mdir.angle==angs(d) & this_mdir.distance==amps(a)), num2cell(this_mdir.SACCADE(this_mdir.angle==angs(d) & this_mdir.distance==amps(a))), 'uni', 1))];
                        end
                        n1_zfr_all = [n1_zfr_all; n1_zfr];
                        n2_zfr_all = [n2_zfr_all; n2_zfr];
                        n12_fr_all = [n12_fr_all; n12_fr];

                        [rho,~] = corr(n1_zfr, n2_zfr);

                        rsc_perDir_mdir(sortIdx(d)) = rho;
                        fr_perDir_mdir(sortIdx(d)) = mean(n12_fr);
                    end

                    n1_fr = cellfun(@(q,w,v) ((sum(q{n1}>w(1) & q{n1}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p1)), this_mdir.TARG_ON, num2cell(this_mdir.SACCADE), 'uni', 1);
                    n2_fr = cellfun(@(q,w,v) ((sum(q{n2}>w(1) & q{n2}<=v))/(v-w(1)))*1000, this_mdir.(sprintf('spiketimes_%d',p2)), this_mdir.TARG_ON, num2cell(this_mdir.SACCADE), 'uni', 1);

                    [rsig_mdir,~] = corr(n1_fr,n2_fr);

                    [rsc_mdir,~] = corr(n1_zfr_all, n1_zfr_all);
                    mnFR_mdir = mean(n12_fr_all);

                    % PURS
                    angs = sort(unique(this_purs.angle));
                    amps = sort(unique(this_purs.pursuitSpeed));

                    if p1_units.pur_sel_dir(n1) > p2_units.pur_sel_dir(n2)
                        bestDir = p1_units.purs_targFR_peakDir(n1);
                    else
                        bestDir = p2_units.purs_targFR_peakDir(n2);
                    end

                    relAngs = mod(angs - bestDir + 180, 360) - 180;
                    [~, sortIdx] = sort(relAngs);

                    [rsc_perDir_purs, fr_perDir_purs] = deal(zeros(1, numel(angs)));
                    n1_zfr_all = []; n2_zfr_all = []; n12_fr_all = [];
                    for d = 1:numel(angs)
                        thisAng = angs(d);

                        n12_fr = []; n1_zfr = []; n2_zfr = [];
                        for a = 1:numel(amps)
                            n12_fr = [n1_fr; [cellfun(@(q,w) ((sum(q{n1}>w & q{n1}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p1))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1); cellfun(@(q,w) ((sum(q{n2}>w & q{n2}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p2))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1)]];
                            
                            n1_zfr = [n1_zfr; zscore(cellfun(@(q,w) ((sum(q{n1}>w & q{n1}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p1))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1))];
                            n2_zfr = [n2_zfr; zscore(cellfun(@(q,w) ((sum(q{n2}>w & q{n2}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p2))(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a)), num2cell(this_purs.PURSUIT_TARG_ON(this_purs.angle==angs(d) & this_purs.pursuitSpeed==amps(a))), 'uni', 1))];
                        end
                        n1_zfr_all = [n1_zfr_all; n1_zfr];
                        n2_zfr_all = [n2_zfr_all; n2_zfr];
                        n12_fr_all = [n12_fr_all; n12_fr];

                        [rho,~] = corr(n1_zfr, n2_zfr);

                        rsc_perDir_purs(sortIdx(d)) = rho;
                        fr_perDir_purs(sortIdx(d)) = mean(n12_fr);
                    end

                    n1_fr = cellfun(@(q,w) ((sum(q{n1}>w & q{n1}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p1)), num2cell(this_purs.PURSUIT_TARG_ON), 'uni', 1);
                    n2_fr = cellfun(@(q,w) ((sum(q{n2}>w & q{n2}<=w+200))/(200))*1000, this_purs.(sprintf('spiketimes_%d',p2)), num2cell(this_purs.PURSUIT_TARG_ON), 'uni', 1);

                    [rsig_purs,~] = corr(n1_fr,n2_fr);

                    [rsc_purs,~] = corr(n1_zfr_all, n1_zfr_all);
                    mnFR_purs = mean(n12_fr_all);
                    
                    pairs(rr,:) = {p1_units.monkey(n1), p1_units.sess_name(n1) ...
                            p1_units.probe_index(n1), p2_units.probe_index(n2), ...
                            p1_units.cluster_id(n1), p2_units.cluster_id(n2), ...
                            p1_units.VMI_quantile(n1), p2_units.VMI_quantile(n2), ...
                            p1_units.SPI_quantile(n1), p2_units.SPI_quantile(n2), ...
                            abs(p1_units.unit_locations{n1}(2)-p2_units.unit_locations{n2}(2)), ...
                            abs(mod(p1_units.vis_pref_dir(n1) - p2_units.vis_pref_dir(n2) + 180, 360) - 180), ...
                            abs(mod(p1_units.sac_pref_dir(n1) - p2_units.sac_pref_dir(n2) + 180, 360) - 180), ...
                            abs(mod(p1_units.pur_pref_dir(n1) - p2_units.pur_pref_dir(n2) + 180, 360) - 180), ...
                            mnFR_mdir, rsc_mdir, rsig_mdir, rsc_perDir_mdir, fr_perDir_mdir, ...
                            mnFR_purs, rsc_purs, rsig_purs, rsc_perDir_purs, fr_perDir_purs
                            };

                    rr = rr + 1;
                end
            end   
            pairs_probes{(p1-1)*2 + p2} = pairs;
        end
    end
    pairs_all{s} = pairs_probes;
end

save('/ix1/pmayo/lab_NHPdata/pairs_table.mat', 'pairs_all');

end