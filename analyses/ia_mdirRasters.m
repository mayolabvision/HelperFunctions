function ia_mdirRasters(data,varargin)
    %UNTITLED2 Summary of this function goes here
    %  8 directions only
    p = inputParser;
    addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'PROBE_INDEX', 1, @isnumeric);
    addParameter(p, 'ALIGN', 'stim', @ischar);
    addParameter(p, 'NET_THRESH', [], @isnumeric);
    addParameter(p, 'X_LIMITS', [-300 500], @isnumeric)
    addParameter(p, 'Y_LIMITS', [], @isnumeric)
    addParameter(p, 'TICK_LENGTH', [], @isnumeric)
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    addParameter(p, 'N_CHUNKS', NaN, @isnumeric);
    addParameter(p, 'CLUSTER', [], @isnumeric);
    addParameter(p, 'IS_FHC', false, @islogical);
    addParameter(p, 'SAVE_PDF', false, @islogical);
    
    parse(p, data, varargin{:});
    data = p.Results.data;
    FIG_PATH = p.Results.FIG_PATH;
    PROBE_INDEX = p.Results.PROBE_INDEX;
    ALIGN = p.Results.ALIGN;
    NET_THRESH = p.Results.NET_THRESH;
    X_LIMITS = p.Results.X_LIMITS;
    Y_LIMITS = p.Results.Y_LIMITS;
    TICK_LENGTH = p.Results.TICK_LENGTH;
    JOB_ID = p.Results.JOB_ID;
    N_CHUNKS = p.Results.N_CHUNKS;
    CLUSTER = p.Results.CLUSTER;
    IS_FHC = p.Results.IS_FHC;
    SAVE_PDF = p.Results.SAVE_PDF;

    fprintf('\n------------------------------\n')
    if ischar(data)
        [~, filename, ~] = fileparts(data);
        load(data,'S');
        fprintf(sprintf('\n----Data loaded for %s----\n',filename))
    else
        S = data;
    end

    clusts_all = S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.cluster_id;
    if ismember('best_channel',S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.Properties.VariableNames)
        chans =  S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.best_channel;
    else
        chans = nan(numel(clusts_all),1);
    end

    if IS_FHC
        units = S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.unit_id;
        scodes = S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.sort_code;
    else
        units = nan(numel(clusts_all),1);
        scodes = nan(numel(clusts_all),1);
    end

    if isempty(CLUSTER)
        if ~isnan(JOB_ID)
            all_units = clusts_all + 1;
            % Split into 50 chunks as a cell array
            chunks = arrayfun(@(i) all_units(...
                floor((i-1)*numel(all_units)/N_CHUNKS)+1 : ...
                floor(i*numel(all_units)/N_CHUNKS)), ...
                1:N_CHUNKS, 'UniformOutput', false);
            ids = (chunks{(JOB_ID+1)});

            clusts = clusts_all(ids);
            chans = chans(ids);
            units = units(ids);
            scodes = scodes(ids);
        else
            clusts = clusts_all;
            chans = chans;
            units = units;
            scodes = scodes;
        end
    else
        clusts = CLUSTER;
        chans = chans(clusts_all==CLUSTER);
        units = units(clusts_all==CLUSTER);
        scodes = scodes(clusts_all==CLUSTER);
    end

    if isempty(chans)
        fprintf('Error: The value given for CLUSTER is higher than the number of clusters recorded.\n');
        fprintf('\n------------------------------\n')
        return;
    end

    probe_label = string(S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.probe_label(1));
    hardware_config = string(S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.hardware_config(1));

    prb_name = sprintf('spiketimes_%d', PROBE_INDEX);

    if ~isempty(FIG_PATH)
        FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config, probe_label), 'mdir_rasters', sprintf('%s_aligned',ALIGN));
        if ~exist(FIG_PATH2, 'dir'), mkdir(FIG_PATH2); end
    else
        FIG_PATH2 = [];
    end

    if isequal(ALIGN,'stim')
        FR_WIN = [50,150];
        xlab = 'time aligned to target onset (ms)';
    elseif isequal(ALIGN,'sacc')
        FR_WIN = [-50,50];
        xlab = 'time aligned to saccade onset (ms)';
    elseif isequal(ALIGN,'fix_off')
        FR_WIN = [0,100];
        xlab = 'time aligned to fixation offset (ms)';

    end

    % Find mdir or dirmem fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'mdir', 'dirmem'}, 'IgnoreCase', true));
    
    T = []; 
    for mm = 1:numel(matchingFields)
        tt = S.(matchingFields{mm}).tbl;
        if isempty(NET_THRESH)
            vars = {'trialName', 'angle', 'distance', 'result', 'TARG_ON', 'SACCADE', 'FIX_OFF', prb_name};
        else
            vars = {'trialName', 'angle', 'distance', 'result', 'TARG_ON', 'SACCADE', 'FIX_OFF', prb_name, sprintf('netlabels_%d', PROBE_INDEX)};
        end
        T = [T; tt(:, vars)];
    end

    T = T(T.result=='CORRECT',:);
    T = T(~cellfun(@(v) sum(cellfun(@(q) isempty(q), v, 'uni', 1)) == numel(T.(prb_name){1}), T.(prb_name), 'uni', 1),:);

    % line_color = {[123,44,191]./255; [230,34,172]./255; [191,44,44]./255};
    % tick_color = {[98,35,152]./255; [184,27,137]./255; [152,35,35]./255};
    % sem_shade = {[228,212,242]./255; [250,210,238]./255; [242,212,212]./255};
    
    line_color = {[42,157,143]./255; [42,114,157]./255; [89,157,42]./255};
    tick_color = {[25,94,85]./255; [25,68,94]./255; [53,94,25]./255};
    sem_shade = {[212,235,232]./255; [212,226,235]./255; [221,235,212]./255};

    for u=1:length(clusts)
        clust = clusts(u);
        unit = units(u);
        scode = scodes(u);

        if IS_FHC
            names = string(T.trialName);
            % Extract the two digits after "unit"
            tokens = regexp(names, 'unit(\d{2})', 'tokens', 'once');
            unitNums = cellfun(@(x) str2double(x), tokens);

            T = T(unitNums == unit, :); 
        end

        T = T(T.result=='CORRECT',:);

        if isempty(T)
            fprintf(sprintf('\n----No mdir trials detected for cluster %d----\n',clust))
        else
             angles = sort(unique(T.angle))';
             if numel(angles)==8
                angle_order = [6,3,2,1,4,7,8,9];
                plotType = '8dir';
             elseif numel(angles)==4
                 angle_order = [6,2,4,8];
                 plotType = '8dir';
             else
                 angle_order = 1:numel(angles);
                 plotType = 'rows';
             end
             distances = sort(unique(T.distance));

            if ~exist(fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, clust, chans(u))), 'file') | isempty(FIG_PATH)
                if ~isempty(FIG_PATH)
                    f3a = figure('Visible','off');
                else
                    f3a = figure('Visible','on');
                end
                f3a.Position = [100 100 1800 900];
            
                y_lims = []; % Store y-axis limits
                frs_perAng = cell(length(angles),length(distances));
                for ang = 1:length(angles)
                    these_trls = T(T.angle==angles(ang),:);
    
                    these_trls.trialName = categorical(regexprep( cellstr(these_trls.trialName), 'mdir1\.(\d+)', 'mdir1.${sprintf(''%04d'', str2double($1))}'));
                    these_trls = sortrows(these_trls, {'distance', 'trialName'});

                    if IS_FHC
                        spks = these_trls.(prb_name);
                        spks = spks(:,scode+1);
                        if isequal(ALIGN,'stim')
                            sptimes = cellfun(@(w,v) w-v(1), spks, these_trls.TARG_ON, 'uni', 0);
                        elseif isequal(ALIGN,'sacc')
                            sptimes = cellfun(@(w,v) w-v(1), spks, num2cell(these_trls.SACCADE), 'uni', 0);
                        elseif isequal(ALIGN,'fix_off')
                            sptimes = cellfun(@(w,v) w-v(1), spks, these_trls.FIX_OFF, 'uni', 0);
                        end

                        if ~isempty(NET_THRESH)
                            nets = these_trls.(sprintf('netlabels_%d', PROBE_INDEX));
                            netlabs = nets(:,scode+1);
                            sptimes = cellfun(@(q,v) q(v>NET_THRESH), sptimes, netlabs, 'uni', 0);
                        end

                    else
                        if isequal(ALIGN,'stim')
                            sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(clust+1)}, these_trls.(prb_name), 'uni', 0), these_trls.TARG_ON, 'uni', 0);
                        elseif isequal(ALIGN,'sacc')
                            sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(clust+1)}, these_trls.(prb_name), 'uni', 0), num2cell(these_trls.SACCADE), 'uni', 0);
                        elseif isequal(ALIGN,'fix_off')
                            sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(clust+1)}, these_trls.(prb_name), 'uni', 0), these_trls.FIX_OFF, 'uni', 0);
                        end

                        if ~isempty(NET_THRESH)
                            netlabs = cellfun(@(q) q{(clust+1)}, these_trls.(sprintf('netlabels_%d', PROBE_INDEX)), 'uni', 0);
                            sptimes = cellfun(@(q,v) q(v>NET_THRESH), sptimes, netlabs, 'uni', 0);
                        end
                    end
    
                    if isequal(plotType,'8dir')
                        subplot(3,3,angle_order(ang))
                    else
                        if numel(angles) > 1
                            subplot(numel(angles)+1,1,angle_order(ang))
                        else
                            subplot(numel(angles),1,angle_order(ang))
                        end
                    end
    
                    if numel(distances)==1
                        if ~isempty(TICK_LENGTH)
                            raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'FR_WINDOW', FR_WIN, 'TICK_LENGTH', TICK_LENGTH)
                        else
                            raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'FR_WINDOW', FR_WIN)
                        end
                        frs_perAng{ang} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/(FR_WIN(2)-FR_WIN(1)))), sptimes, 'uni', 1);
                    else
                        [line_colors, tick_colors, sem_shades] = deal(cell(height(these_trls),1)); 
                        for dd = 1:numel(distances)
                            line_colors(these_trls.distance==distances(dd)) = line_color(dd);
                            tick_colors(these_trls.distance==distances(dd)) = tick_color(dd);
                            sem_shades(these_trls.distance==distances(dd)) = sem_shade(dd);
    
                            frs_perAng{ang,dd} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/(FR_WIN(2)-FR_WIN(1)))), sptimes(these_trls.distance==distances(dd)), 'uni', 1);
                        end
    
                        if ~isempty(TICK_LENGTH)
                            raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_colors, 'TICK_COLOR', tick_colors, 'FR_WINDOW', FR_WIN, 'TICK_LENGTH', TICK_LENGTH)
                        else
                            raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_colors, 'TICK_COLOR', tick_colors, 'FR_WINDOW', FR_WIN)
                        end
                    end
            
                    yyaxis left;
                    ax = gca;
                    y_lims = [y_lims; ax.YLim];
            
                end
            
                % Pad each cell with NaNs to match maxLength
                maxLength = max(cellfun(@numel, frs_perAng)); maxLength = max(maxLength);
    
                str_title = deal(cell(1,numel(distances)));
                for dd = 1:length(distances)
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
                    theta = angles;
                    [visds, visdp] = tuningbias(theta,mean(stimrate,'omitnan'));
                    
                    if numel(angles) > 1
                        if dd==1
                            if isequal(plotType,'8dir')
                                subplot(3,3,5)
                            else
                                subplot(numel(angles)+1,1,numel(angles)+1)
                            end
                        end
                        rho = mean(stimrate,'omitnan');
                        h1(dd) = polarplot(deg2rad([theta 0]),[rho rho(1)],'o-',...
                            'markerfacecolor',line_color{dd},'linewidth',3,'color',line_color{dd});
                        hold on
                        polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'o--','LineWidth',2,'Color',sem_shade{dd});
                        polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'o--','LineWidth',2,'Color',sem_shade{dd});   
                        
                        polarplot(deg2rad(visdp), max(rho), '^', 'MarkerFaceColor', line_color{dd}, 'MarkerEdgeColor', line_color{dd}, 'MarkerSize', 10);
        
                        tcolor = tick_color{dd};
                        str_title{dd} = sprintf('%d deg -- Dir: %0.2f, Sel: %0.2f', distances(dd), visdp, visds);
                    end

                    if dd==length(distances)
                         title(str_title);
                    end
                end 
   
                legend_labels = arrayfun(@(d) sprintf('%d deg', distances(d)), 1:length(distances), 'UniformOutput', false);
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
                    if isequal(plotType,'8dir')
                        subplot(3,3,ang)
                        if ang ~= 5
                            yyaxis left;
                            ylim(global_y_lim);
                        end
                    else
                        if numel(angles) > 1
                            subplot(numel(angles),1,angle_order(ang))
                            yyaxis left;
                            ylim(global_y_lim);
                        end
                    end
                end
    
                han=axes(f3a,'visible','off'); 
                han.Title.Visible='on';
                han.XLabel.Visible='on';
                xlabel(han,{'';xlab},'fontsize',16);
    
                title(han, {
                    sprintf('%s --- %s --- cluster %d (channel %d)',S.sess_name, probe_label, clust, chans(u));
                    '' 
                }, 'fontsize',16,'interpreter','none')
            
                if ~isempty(FIG_PATH)
                    if SAVE_PDF
                        savebigPDF(f3a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, clust, chans(u))));
                    else
                        print(f3a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, clust, chans(u))), '-dpng', '-r200');
                    end
                end
    
                fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE----',PROBE_INDEX, clust))
    
            else
                fprintf(sprintf('\n----PROBE %d, Unit %.4d exists----',PROBE_INDEX, clust))
            end
        end
    end
    fprintf('\n------------------------------\n')
end
