function ia_pursRasters(data,varargin)
    % ff
    % ff

    p = inputParser;
    addRequired(p, 'data', @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'PROBE_INDEX', 1, @isnumeric);
    addParameter(p, 'ALIGN', 'targ', @ischar);
    addParameter(p, 'PURE_ONLY', false, @islogical)
    addParameter(p, 'X_LIMITS', [-300 500], @isnumeric)
    addParameter(p, 'Y_LIMITS', [], @isnumeric)
    addParameter(p, 'TICK_LENGTH', [], @isnumeric)
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    addParameter(p, 'N_CHUNKS', NaN, @isnumeric);
    addParameter(p, 'CLUSTER', [], @isnumeric);
    addParameter(p, 'SAVE_PDF', false, @islogical);   
 
    parse(p, data, varargin{:});
    data = p.Results.data;
    FIG_PATH = p.Results.FIG_PATH;
    PROBE_INDEX = p.Results.PROBE_INDEX;
    ALIGN = p.Results.ALIGN;
    PURE_ONLY = p.Results.PURE_ONLY;
    X_LIMITS = p.Results.X_LIMITS;
    Y_LIMITS = p.Results.Y_LIMITS;
    TICK_LENGTH = p.Results.TICK_LENGTH;
    JOB_ID = p.Results.JOB_ID;
    N_CHUNKS = p.Results.N_CHUNKS;
    CLUSTER = p.Results.CLUSTER;
    SAVE_PDF = p.Results.SAVE_PDF;

    fprintf('\n------------------------------\n')
    if ischar(data)
        %[~, filename, ~] = fileparts(data);
        load(data, 'S');
        fprintf('\n----Data loaded----\n')
    else
        %filename = data.sessionName;
        S = data;
    end

    units = S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.cluster_id;
    if ismember('best_channel',S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.Properties.VariableNames)
        chans =  S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.best_channel;
    else
        chans = nan(numel(units),1);
    end

    if isempty(CLUSTER)
        snrs  = S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.snr;
        depths = cellfun(@(q) q(2), S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.unit_locations, 'uni', 1)./1000;
        kslabs = S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.KSLabel_clusters;

        if ~isnan(JOB_ID)
            all_units = units + 1;
            % Split into 50 chunks as a cell array
            chunks = arrayfun(@(i) all_units(...
                floor((i-1)*numel(all_units)/N_CHUNKS)+1 : ...
                floor(i*numel(all_units)/N_CHUNKS)), ...
                1:N_CHUNKS, 'UniformOutput', false);
            ids = (chunks{(JOB_ID+1)});

            units = units(ids);
            chans = chans(ids);

            snrs = snrs(ids);
            depths = depths(ids);
            kslabs = kslabs(ids);

        end
    else
        units = CLUSTER;
        chans = chans(units==CLUSTER);
    end

    probe_label = string(S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.probe_label(1));
    hardware_config = string(S.kilosort([S.kilosort.probe_index] == PROBE_INDEX).clusters.hardware_config(1));

    if ~isempty(FIG_PATH)
        if PURE_ONLY
            FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config, probe_label), 'purs_rasters', sprintf('pure_only-%s_aligned',ALIGN));
        else
            FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config, probe_label), 'purs_rasters', sprintf('all_trls-%s_aligned',ALIGN));
        end
        if ~exist(FIG_PATH2, 'dir'), mkdir(FIG_PATH2); end
    else
        FIG_PATH2 = 'blahblah';
    end

    if isequal(ALIGN,'targ')
        FR_WIN = [50, 250];
        xlab = 'time aligned to target motion onset (ms)';
    elseif isequal(ALIGN,'purs')
        FR_WIN = [-50,50];
        xlab = 'time aligned to pursuit onset (ms)';
    end

    % Find rfmp or rfMapping fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'purs', 'pursuit'}, 'IgnoreCase', true));
    
    if ~isempty(FIG_PATH)
        if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end   
    end

    T = []; 
    for mm = 1:numel(matchingFields)
        T = [T; S.(matchingFields{mm}).tbl];
    end

    T = T(T.result=='CORRECT',:);

    if ~ismember('eyePos',T.Properties.VariableNames)
         eyePos = cellfun(@(x) filterEyeTraces_EyeLink(x), T.eyedata, 'uni', 0);
        [eyeVel, eyeAcc] = cellfun(@(x) calcDerivative_eyeTraces(x), eyePos, 'uni', 0);
    
        T.eyePos = eyePos; T.eyeVel = eyeVel; T.eyeAcc = eyeAcc;
    end

    if ~ismember('saccades', T.Properties.VariableNames)
        saccades = cellfun(@(q) detect_saccades(q, 'VEL_THRESH', 30, 'ACC_THRESH', 500), cellfun(@(v,f) v(:,f(1):end), eyeVel, num2cell(T.FIXATE), 'uni', 0), 'uni', 0);
        T.saccades = cellfun(@(q,v) num2cell(q+v(1),2), saccades, num2cell(T.FIXATE), 'uni', 0);
    end

    if ~ismember('pursuitOnset', T.Properties.VariableNames)
        [pursuitOnset, pursuitLatency] = cellfun(@(u,v,w) detect_pursuitOnset(u, v, w, 'PLOT_TRACES', false), T.eyeVel, num2cell(T.PURSUIT_TARG_ON), num2cell(T.pursuitSpeed), 'uni', 1); 
        T.pursuitOnset = pursuitOnset;
        T.pursuitLatency = pursuitLatency;
    end

    if ~ismember('pursType', T.Properties.VariableNames)
        csTrials = cellfun(@(u,v) sum(((cellfun(@(q) q(2), u, 'uni', 1) - (v)) >= 0) & ((cellfun(@(q) q(1), u, 'uni', 1) - (v+T.params(1).block.crossingTime)) < 150)), T.saccades, num2cell(T.PURSUIT_TARG_ON), 'uni', 1);
        T.pursType = repmat(categorical("sacc"),height(T),1);
        T.pursType(~logical(csTrials)) = "pure";
    end

    %T = T(T.msOffset<0 | isnan(T.msOffset),:);

    if PURE_ONLY
        T = T(T.pursType=='pure' & T.jump==-1,:);
    end
    
    angles = sort(unique(T.angle))';
    angle_order = [6,3,2,1,4,7,8,9];
    speeds = sort(unique(T.pursuitSpeed));

    if PROBE_INDEX==1 % purples/pinks
        line_color = {[123,44,191]./255; [230,34,172]./255};
        tick_color = {[98,35,152]./255; [184,27,137]./255};
        sem_shade = {[228,212,242]./255; [250,210,238]./255};
    else % greens/blues
        line_color = {[42,157,143]./255; [42,114,157]./255};
        tick_color = {[25,94,85]./255; [25,68,94]./255};
        sem_shade = {[212,235,232]./255; [212,226,235]./255};
    end

    prb_name = ['spiketimes_' num2str(PROBE_INDEX)];
    for u=1:length(units)
        unit = units(u);
        if ~exist(fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), 'file') | isempty(FIG_PATH)
            if ~isempty(FIG_PATH)
                f3a = figure('Visible','off');
            else
                f3a = figure('Visible','on');
            end
            f3a.Position = [100 100 1800 900];
        
            y_lims = []; % Store y-axis limits
            frs_perAng = cell(length(angles),length(speeds));
            for ang = 1:length(angles)
                these_trls = T(T.angle==angles(ang),:);

                these_trls.trialName = categorical(regexprep( cellstr(these_trls.trialName), 'purs1\.(\d+)', 'purs1.${sprintf(''%04d'', str2double($1))}'));
                these_trls = sortrows(these_trls, {'pursuitSpeed', 'trialName'});

                if isequal(ALIGN,'targ')
                    sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(prb_name), 'uni', 0), num2cell(these_trls.PURSUIT_TARG_ON), 'uni', 0);
                elseif isequal(ALIGN,'purs')
                    sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(prb_name), 'uni', 0), num2cell(these_trls.pursuitOnset), 'uni', 0);
                end

                subplot(3,3,angle_order(ang))

                if numel(speeds)==1
                    if ~isempty(TICK_LENGTH)
                        raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'SEM_SHADE', sem_shade{1}, 'FR_WINDOW', FR_WIN, 'TICK_LENGTH', TICK_LENGTH)
                    else
                        raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'SEM_SHADE', sem_shade{1}, 'FR_WINDOW', FR_WIN)
                    end
                    frs_perAng{ang} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/(FR_WIN(2)-FR_WIN(1)))), sptimes, 'uni', 1);
                else
                    [line_colors, tick_colors, sem_shades] = deal(cell(height(these_trls),1)); 
                    for dd = 1:numel(speeds)
                        line_colors(these_trls.pursuitSpeed==speeds(dd)) = line_color(dd);
                        tick_colors(these_trls.pursuitSpeed==speeds(dd)) = tick_color(dd);
                        sem_shades(these_trls.pursuitSpeed==speeds(dd)) = sem_shade(dd);

                        frs_perAng{ang,dd} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/(FR_WIN(2)-FR_WIN(1)))), sptimes(these_trls.pursuitSpeed==speeds(dd)), 'uni', 1);
                    end

                    if size(sptimes,1) > 1
                        if ~isempty(TICK_LENGTH)
                            raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_colors, 'SEM_SHADE', sem_shades, 'TICK_COLOR', tick_colors, 'FR_WINDOW', FR_WIN, 'TICK_LENGTH', TICK_LENGTH)
                        else
                            raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_colors, 'SEM_SHADE', sem_shades, 'TICK_COLOR', tick_colors, 'FR_WINDOW', FR_WIN)
                        end
                    end
                end
        
                yyaxis left;
                ax = gca;
                y_lims = [y_lims; ax.YLim];
        
            end
        
            % Pad each cell with NaNs to match maxLength
            maxLength = max(cellfun(@numel, frs_perAng)); maxLength = max(maxLength);

            str_title = deal(cell(1,numel(speeds)));
            for dd = 1:length(speeds)
                frs_perAng2 = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_perAng(:,dd), 'UniformOutput', false);
                
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
                %theta = 0:360/length(a.CND):360; theta(end)=[];
                theta = 0:45:315;
                [visds, visdp] = tuningbias(theta,mean(stimrate,'omitnan'));
                
                if dd==1
                    subplot(3,3,5)
                end
                rho = mean(stimrate,'omitnan');
                h1(dd) = polarplot(deg2rad([theta 0]),[rho rho(1)],'o-',...
                    'markerfacecolor',line_color{dd},'linewidth',3,'color',line_color{dd});
                hold on
                polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'o--','LineWidth',2,'Color',sem_shade{dd});
                polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'o--','LineWidth',2,'Color',sem_shade{dd});   
                
                polarplot(deg2rad(visdp), max(rho), '^', 'MarkerFaceColor', line_color{dd}, 'MarkerEdgeColor', line_color{dd}, 'MarkerSize', 10);

                tcolor = tick_color{dd};
                str_title{dd} = sprintf('%d deg/s -- Dir: %0.2f, Sel: %0.2f', speeds(dd), visdp, visds);
            end 

            title(str_title);
               
            legend_labels = arrayfun(@(d) sprintf('%d deg/s', speeds(d)), 1:length(speeds), 'UniformOutput', false);
            legend(h1, legend_labels, 'Location', 'best');
            prettyFig;


            % Find the global y-axis limits
            if isempty(Y_LIMITS)
                global_y_lim = [min(y_lims(:,1)), max(y_lims(:,2))];
            else
                global_y_lim = Y_LIMITS;
            end
            
            % Apply the limits to all subplots
            for ang = 1:max(angle_order)
                subplot(3,3,ang)
                if ang ~= 5
                    yyaxis left;
                    ylim(global_y_lim);
                end
            end

            han=axes(f3a,'visible','off'); 
            han.Title.Visible='on';
            han.XLabel.Visible='on';
            xlabel(han,{'';xlab},'fontsize',16);

            title(han, {
                sprintf('%s --- %s --- cluster %d (channel %d)',S.sess_name, probe_label, unit, chans(u));
                sprintf('ks_label = %s, snr = %.4f, y_pos = %.2f um', kslabs{u}, snrs(u), depths(u)) 
            }, 'fontsize',16,'interpreter','none')

            if ~isempty(FIG_PATH)
                if SAVE_PDF
                    savebigPDF(f3a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))));
                else
                    print(f3a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), '-dpng', '-r200');
                end
            end
            fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE----',PROBE_INDEX, unit))

        else
            fprintf(sprintf('\n----PROBE %d, Unit %.4d exists----',PROBE_INDEX, unit))
        end

    end
    fprintf('\n------------------------------\n')
end
