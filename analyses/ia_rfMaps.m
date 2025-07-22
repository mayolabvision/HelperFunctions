function ia_rfMaps(data,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'IMEC', 0, @isnumeric);
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    addParameter(p, 'N_CHUNKS', NaN, @isnumeric);
    addParameter(p, 'CLUSTER', [], @isnumeric);
    
    parse(p, data, varargin{:});
    data = p.Results.data;
    FIG_PATH = p.Results.FIG_PATH;
    IMEC = p.Results.IMEC;
    JOB_ID = p.Results.JOB_ID;
    N_CHUNKS = p.Results.N_CHUNKS;
    CLUSTER = p.Results.CLUSTER;

    fprintf('\n------------------------------\n')
    if ischar(data)
        [~, filename, ~] = fileparts(data);
        load(data,'S');
        fprintf(sprintf('\n----Data loaded for %s----\n',filename))
    else
        filename = data.sess_name;
        S = data;
    end

    units = S.kilosort([S.kilosort.imec] == IMEC).clusters.cluster_id; 
    if ismember('channel_id',S.kilosort([S.kilosort.imec] == IMEC).clusters.Properties.VariableNames)
        chans =  S.kilosort([S.kilosort.imec] == IMEC).clusters.channel_id;
    else
        chans = nan(numel(units),1);
    end

    if isempty(CLUSTER)
        % snrs  = S.kilosort([S.kilosort.imec] == IMEC).clusters.snr;
        % depths = S.kilosort([S.kilosort.imec] == IMEC).clusters.y_pos;
        % kslabs = S.kilosort([S.kilosort.imec] == IMEC).clusters.KSLabel_clusters;

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
         
            % snrs = snrs(ids);
            % depths = depths(ids);
            % kslabs = kslabs(ids);
            
        end
    else
        units = CLUSTER;
        chans = chans(units==CLUSTER);
    end
    
    % Find rfmp or rfMapping fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'rfmp', 'rfMapping'}, 'IgnoreCase', true));

    task_conditions = zeros(length(matchingFields),2);
    for mm = 1:numel(matchingFields)
        stimuli = cell2mat(vertcat(S.(matchingFields{mm}).tbl.conditions{:}));
        task_conditions(mm,:) = [length(sort(unique(stimuli(:,1))).') * length(sort(unique(stimuli(:,2))).') S.(matchingFields{mm}).params.frameCount(1)];
    end

    [unique_conditions, ~, row_conds] = unique(task_conditions, 'rows', 'stable');

    if numel(unique(row_conds))==1
        if ~isempty(FIG_PATH)
            if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end    
        end
        T = []; 
        for mm = 1:numel(matchingFields)
            T = [T; S.(matchingFields{mm}).tbl];
        end

        for u=1:length(units)
            unit = units(u); 
            if ~exist(fullfile(FIG_PATH, sprintf('imec%d_unit%04d_chan%03d.png', IMEC, unit, chans(u))), 'file')  | isempty(FIG_PATH)
                [frs,bin_edges,xvals,yvals] = format_tableToRFMap(T, 'IMEC', IMEC, 'UNITS', (unit+1));
            
                f2a = figure('Visible','off');
                f2a.Position = [100 100 1800 900];
                tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',xvals, 'Y_VALS',yvals);
                
                if (IMEC)==0
                    title(tl,sprintf('%s --- LEFT --- cluster %d (channel %d)',filename, unit, chans(u)),'fontsize',16,'interpreter','none')
                else
                    title(tl,sprintf('%s --- RIGHT --- cluster %d (channel %d)',filename, unit, chans(u)),'fontsize',16,'interpreter','none')
                end
                %subtitle(tl, sprintf('ks_label = %s, snr = %.4f, y_pos = %.2f um', kslabs(u), snrs(u), depths(u)),'fontsize',12,'interpreter','none')
                
                annotation('textbox', [0.8 0.89 0.2 0.1], ... % [x y w h] in normalized figure units
                    'String', sprintf('N = %d repeats', min(min(min(cellfun(@length, frs{1}))))), ...
                    'FontSize', 16, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'right');
                
                if ~isempty(FIG_PATH)
                    print(f2a, fullfile(FIG_PATH, sprintf('imec%d_unit%04d_chan%03d.png', IMEC, unit, chans(u))), '-dpng', '-r200');
                end
                fprintf(sprintf('\n----IMEC %d, Unit %.4d COMPLETE----',IMEC, unit))
            else
                fprintf(sprintf('\n----IMEC %d, Unit %.4d exists----',IMEC, unit))
            end
        end
    else
        for mm = 1:size(unique_conditions,1)
            these_rows = find(row_conds==mm);
            T = []; 
            for r = 1:numel(these_rows)
                T = [T; S.(matchingFields{these_rows(r)}).tbl];
            end

            stim_duration_ms = T.params(1,1).block.frameCount*(1000/120);
            
            if ~isempty(FIG_PATH)
                if ~exist(fullfile(FIG_PATH, strjoin(matchingFields(these_rows), '_')), 'dir'), mkdir(fullfile(FIG_PATH, strjoin(matchingFields(these_rows), '_'))); end    
            end

            for u=1:length(units)
                unit = units(u); 
                if ~exist(fullfile(FIG_PATH, strjoin(matchingFields(these_rows), '_'), sprintf('imec%d_unit%04d_chan%03d.png', IMEC, unit, chans(u))), 'file') | ~isempty(FIG_PATH)
                    [frs,bin_edges,xvals,yvals] = format_tableToRFMap(T, 'IMEC', IMEC, 'UNITS', (unit+1));
                
                    f2a = figure('Visible','off');
                    f2a.Position = [100 100 1800 900];
                    tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',xvals, 'Y_VALS',yvals);
                    
                    if (IMEC)==0
                        title(tl,sprintf('%s_%s --- LEFT --- cluster %d (channel %d)',filename, strjoin(matchingFields(these_rows), '_'), unit, chans(u)),'fontsize',16,'interpreter','none')
                    else
                        title(tl,sprintf('%s_%s --- RIGHT --- cluster %d (channel %d)',filename, strjoin(matchingFields(these_rows), '_'), unit, chans(u)),'fontsize',16,'interpreter','none')
                    end
                    %subtitle(tl, sprintf('ks_label = %s, snr = %.4f, y_pos = %.2f um', kslabs(u), snrs(u), depths(u)),'fontsize',12,'interpreter','none')
                    
                    annotation('textbox', [0.77 0.89 0.2 0.1], ... % [x y w h] in normalized figure units
                        'String', sprintf('N = %d repeats', min(min(min(cellfun(@length, frs{1}))))), ...
                        'FontSize', 16, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'right');

                    annotation('textbox', [0.01 0.89 0.2 0.1], ... % [x y w h] in normalized figure units
                        'String', sprintf('probe duration = %.1f ms', stim_duration_ms), ...
                        'FontSize', 16, ...
                        'EdgeColor', 'none', ...
                        'HorizontalAlignment', 'left');
                    

                    if ~isempty(FIG_PATH)
                        print(f2a, fullfile(FIG_PATH, strjoin(matchingFields(these_rows), '_'), sprintf('imec%d_unit%04d_chan%03d.png', IMEC, unit, chans(u))), '-dpng', '-r200');
                    end
                    fprintf(sprintf('\n----IMEC %d, Unit %.4d COMPLETE----',IMEC, unit))
                else
                    fprintf(sprintf('\n----IMEC %d, Unit %.4d exists----',IMEC, unit))
                end
            end
        end
    end
    fprintf('\n------------------------------\n')
end
