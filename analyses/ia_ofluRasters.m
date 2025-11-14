function ia_ofluRasters(data,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'PROBE_INDEX', 1, @isnumeric);
    addParameter(p, 'ALIGN', 'stim', @ischar);
    addParameter(p, 'X_LIMITS', [-100 700], @isnumeric)
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
    X_LIMITS = p.Results.X_LIMITS;
    Y_LIMITS = p.Results.Y_LIMITS;
    TICK_LENGTH = p.Results.TICK_LENGTH;
    JOB_ID = p.Results.JOB_ID;
    N_CHUNKS = p.Results.N_CHUNKS;
    CLUSTER = p.Results.CLUSTER;
    SAVE_PDF = p.Results.SAVE_PDF;

    fprintf('\n------------------------------\n')
    if ischar(data)
        [~, filename, ~] = fileparts(data);
        load(data,'S');
        fprintf(sprintf('\n----Data loaded for %s----\n',filename))
    else
        S = data;
    end

    units = S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.cluster_id;
    if ismember('best_channel',S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.Properties.VariableNames)
        chans =  S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.best_channel;
    else
        chans = nan(numel(units),1);
    end

    if isempty(CLUSTER)

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
        end
    else
        units = CLUSTER;
        chans = chans(units==CLUSTER);
    end

    probe_label = string(S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.probe_label(1));
    hardware_config = string(S.sorting([S.sorting.probe_index] == PROBE_INDEX).clusters.hardware_config(1));

    prb_name = sprintf('spiketimes_%d',PROBE_INDEX);

    if ~isempty(FIG_PATH)
        FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config, probe_label), 'oflu_rasters', sprintf('%s_aligned', ALIGN));
        if ~exist(FIG_PATH2, 'dir'), mkdir(FIG_PATH2); end
    else
        FIG_PATH2 = [];
    end

    if isequal(ALIGN,'stim')
        FR_WIN = [0,515];
        xlab = 'time aligned to stim onset (ms)';
    end

    % Find mdir or dirmem fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'oflu'}, 'IgnoreCase', true));
    
    T = []; 
    for mm = 1:numel(matchingFields)
        tt = S.(matchingFields{mm}).tbl;
        vars = {'recColor', 'result', 'STIM_ON', 'STIM_OFF', 'trialName', prb_name};
        T = [T; tt(:, vars)];
    end

    T = T(T.result=='CORRECT',:);
    
    recColors = sort(unique(T.recColor))';

    if PROBE_INDEX==1 % purples/pinks
        line_color = {[123,44,191]./255; [230,34,172]./255; [191,44,44]./255};
        tick_color = {[98,35,152]./255; [184,27,137]./255; [152,35,35]./255};
        sem_shade = {[228,212,242]./255; [250,210,238]./255; [242,212,212]./255};
    else % greens/blues
        line_color = {[42,157,143]./255; [42,114,157]./255; [89,157,42]./255};
        tick_color = {[25,94,85]./255; [25,68,94]./255}; [53,94,25]./255;
        sem_shade = {[212,235,232]./255; [212,226,235]./255; [221,235,212]./255};
    end

    for u=1:length(units)
        unit = units(u);
        if ~exist(fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), 'file') | isempty(FIG_PATH)
            if ~isempty(FIG_PATH)
                f3a = figure; %('Visible','off');
            else
                f3a = figure; %('Visible','on');
            end
            f3a.Position = [100 100 1100 900];
            
            y_lims = [];
            for col = 1:length(recColors)
                subplot(numel(recColors),1,col)
                these_trls = T(T.recColor==recColors(col),:);

                these_trls.trialName = categorical(regexprep( cellstr(these_trls.trialName), 'oflu\.(\d+)', 'oflu.${sprintf(''%04d'', str2double($1))}'));
                these_trls = sortrows(these_trls, {'recColor', 'trialName'});

                if isequal(ALIGN,'stim')
                    if size(S.sorting,1) == 1
                        sptimes = cellfun(@(w,v) w-v(1), these_trls.(prb_name), num2cell(these_trls.STIM_ON), 'uni', 0);
                    else
                        sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(prb_name), 'uni', 0), num2cell(these_trls.STIM_ON), 'uni', 0);
                    end
                end

                netlabs = cellfun(@(q) q{(unit+1)}, these_trls.net_labels, 'uni', 0);
                sptimes = cellfun(@(q,v) q(v>0.02), sptimes, netlabs, 'uni', 0);

                if ~isempty(TICK_LENGTH)
                    raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'SEM_SHADE', sem_shade{1}, 'FR_WINDOW', FR_WIN, 'TICK_LENGTH', TICK_LENGTH)
                else
                    raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'SEM_SHADE', sem_shade{1}, 'FR_WINDOW', FR_WIN)
                end   
        
                yyaxis left;
                ax = gca;
                y_lims = [y_lims; ax.YLim];

                title(sprintf('[128 128 128] --> [%d %d %d]', recColors(col), recColors(col), recColors(col)))
            end
        
            prettyFig;

            % Find the global y-axis limits
            if isempty(Y_LIMITS)
                global_y_lim = [min(y_lims(:,1)), max(y_lims(:,2))];
            else
                global_y_lim = Y_LIMITS;
            end
            
            % Apply the limits to all subplots
            for col = 1:numel(recColors)
                subplot(numel(recColors),1,col)
                yyaxis left;
                ylim(global_y_lim);
            end

            han=axes(f3a,'visible','off'); 
            han.Title.Visible='on';
            han.XLabel.Visible='on';
            xlabel(han,{'';xlab},'fontsize',16);

            title(han, {
            sprintf('%s --- %s --- unit %d (channel %d)',S.sess_name, probe_label, unit, chans(u));
            ''
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
