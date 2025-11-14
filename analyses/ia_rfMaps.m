function ia_rfMaps(data,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data',  @(x) (ischar(x)) || isstruct(x));
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'PROBE_INDEX', 1, @isnumeric);
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    addParameter(p, 'N_CHUNKS', NaN, @isnumeric);
    addParameter(p, 'CLUSTER', [], @isnumeric);
    addParameter(p, 'SAVE_PDF', false, @islogical);
    
    parse(p, data, varargin{:});
    data = p.Results.data;
    FIG_PATH = p.Results.FIG_PATH;
    PROBE_INDEX = p.Results.PROBE_INDEX;
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
        filename = data.sess_name;
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
    
    % Find rfmp or rfMapping fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'rfmp', 'rfMapping'}, 'IgnoreCase', true));

    matchingFields = matchingFields(cellfun(@(f) isfield(S.(f), 'tbl') && ~isempty(S.(f).tbl), matchingFields));

    task_conditions = cell(length(matchingFields), 5);
    for mm = 1:numel(matchingFields)
        stimuli = cell2mat(vertcat(S.(matchingFields{mm}).tbl.conditions{:}));
        task_conditions{mm, 1} = length(unique(stimuli(:,1))) * length(unique(stimuli(:,2))); 
        task_conditions{mm, 2} = S.(matchingFields{mm}).params.frameCount(1); 
        task_conditions{mm, 3} = S.(matchingFields{mm}).params.nStimPerFix(1);
        task_conditions{mm, 4} = S.(matchingFields{mm}).params.bgColor(1,:);
        task_conditions{mm, 5} = [S.(matchingFields{mm}).params.colorR S.(matchingFields{mm}).params.colorG S.(matchingFields{mm}).params.colorB];
    end

    row_strings = cell(size(task_conditions,1),1);
    for i = 1:size(task_conditions,1)
        row_strings{i} = strjoin(cellfun(@mat2str, task_conditions(i,:), 'UniformOutput', false), '|');
    end
    [~, unique_idx, row_conds] = unique(row_strings, 'stable');
    unique_conditions = task_conditions(unique_idx, :);

    if numel(unique(row_conds))==1
        if ~isempty(FIG_PATH)
            FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config, probe_label), 'rfmp_heatmaps'); 
            if ~exist(FIG_PATH2, 'dir'), mkdir(FIG_PATH2); end    
        end
        T = []; 
        for mm = 1:numel(matchingFields)
            T = [T; S.(matchingFields{mm}).tbl];
        end

        T.STIM_ON(T.result~='CORRECT') = cellfun(@(q) q(1:end-1), T.STIM_ON(T.result~='CORRECT'), 'uni', 0);
        T.conditions(T.result~="CORRECT") = cellfun(@(q) q(1:end-1), T.conditions(T.result~='CORRECT'), 'uni', 0);
        T = T(~cellfun(@(q) any(isnan(q)), T.STIM_OFF, 'uni', 1),:);

        stim_duration_ms = T.params(1,1).block.frameCount*(1000/120);
        probe_rgb = [T.params(1,1).block.colorR, T.params(1,1).block.colorG, T.params(1,1).block.colorB];
        bg_rgb = T.params(1,1).block.bgColor;

        bg_text = '';
        if isequal(bg_rgb, [0 0 0])
            bg_text = 'black bg';
        elseif isequal(bg_rgb, [255 255 255])
            bg_text = 'white bg';
        end

        for u=1:length(units)
            unit = units(u); 
            if ~exist(fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), 'file')  | isempty(FIG_PATH)
                [frs,bin_edges,xvals,yvals] = format_tableToRFMap(T, 'PROBE_INDEX', PROBE_INDEX, 'UNITS', (unit+1));
            
                f2a = figure('Visible','off');
                f2a.Position = [100 100 1800 900];
                tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',xvals, 'Y_VALS',yvals,'PROBE_DUR',stim_duration_ms);
                

                title(tl,sprintf('%s --- %s --- cluster %d (channel %d)',S.sess_name, probe_label, unit, chans(u)),'fontsize',16,'interpreter','none')
                subtitle(tl, sprintf('ks_label = %s, snr = %.4f, y_pos = %.2f um', kslabs{u}, snrs(u), depths(u)),'fontsize',12,'interpreter','none')
                
                annotation('textbox', [0.77 0.89 0.2 0.1], ... % [x y w h] in normalized figure units
                           'String', sprintf('N = %d repeats\n%s', min(min(min(cellfun(@length, frs{1})))), bg_text), ...
                           'FontSize', 14, ...
                           'EdgeColor', 'none', ...
                           'HorizontalAlignment', 'right');
                
                if ~isempty(FIG_PATH)
                    print(f2a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), '-dpng', '-r200');
                end
                fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE----',PROBE_INDEX, unit))
            else
                fprintf(sprintf('\n----PROBE %d, Unit %.4d exists----',PROBE_INDEX, unit))
            end
        end
    else
        for mm = 1:size(unique_conditions,1)
            these_rows = find(row_conds==mm);
            T = []; 
            for r = 1:numel(these_rows)
                T = [T; S.(matchingFields{these_rows(r)}).tbl];
            end
            
            T.STIM_ON(T.result~='CORRECT') = cellfun(@(q) q(1:end-1), T.STIM_ON(T.result~='CORRECT'), 'uni', 0);
            T.conditions(T.result~="CORRECT") = cellfun(@(q) q(1:end-1), T.conditions(T.result~='CORRECT'), 'uni', 0);
            T = T(~cellfun(@(q) any(isnan(q)), T.STIM_OFF, 'uni', 1),:);

            stim_duration_ms = T.params(1,1).block.frameCount*(1000/120);
            probe_rgb = [T.params(1,1).block.colorR, T.params(1,1).block.colorG, T.params(1,1).block.colorB];
            bg_rgb = T.params(1,1).block.bgColor;

            bg_text = '';
            if isequal(bg_rgb, [0 0 0])
                bg_text = 'black bg';
            elseif isequal(bg_rgb, [255 255 255])
                bg_text = 'white bg';
            end

            if ~isempty(FIG_PATH)
                FIG_PATH2 = fullfile(FIG_PATH, sprintf('%s_%s',hardware_config, probe_label), 'rfmp_heatmaps', strjoin(matchingFields(these_rows), '_')); 
                if ~exist(FIG_PATH2, 'dir'), mkdir(FIG_PATH2); end    
            end

            for u=1:length(units)
                unit = units(u); 
                if ~exist(fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), 'file') | ~isempty(FIG_PATH)
                    [frs,bin_edges,xvals,yvals] = format_tableToRFMap(T, 'PROBE_INDEX', PROBE_INDEX, 'UNITS', (unit+1));

                    f2a = figure('Visible','off');
                    f2a.Position = [100 100 1800 900];
                    tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',xvals, 'Y_VALS',yvals,'PROBE_DUR',stim_duration_ms);
                    
                    title(tl,sprintf('%s --- %s --- %s --- cluster %d (channel %d)',S.sess_name, probe_label, strjoin(matchingFields(these_rows), '_'), unit, chans(u)),'fontsize',16,'interpreter','none')
                    subtitle(tl, sprintf('ks_label = %s, snr = %.4f, y_pos = %.2f um', kslabs{u}, snrs(u), depths(u)),'fontsize',12,'interpreter','none')
                  
                    
                    annotation('textbox', [0.77 0.89 0.2 0.1], ... % [x y w h] in normalized figure units
                               'String', sprintf('N = %d repeats\n%s', min(min(min(cellfun(@length, frs{1})))), bg_text), ...
                               'FontSize', 14, ...
                               'EdgeColor', 'none', ...
                               'HorizontalAlignment', 'right');

                    
                    if ~isempty(FIG_PATH)
                        if SAVE_PDF
                            savebigPDF(f2a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))));
                        else
                            print(f2a, fullfile(FIG_PATH2, sprintf('%s_clust%04d_chan%03d.png', probe_label, unit, chans(u))), '-dpng', '-r200');
                        end
                    end
                    fprintf(sprintf('\n----PROBE %d, Unit %.4d COMPLETE----',PROBE_INDEX, unit))
                else
                    fprintf(sprintf('\n----PROBE %d, Unit %.4d exists----',PROBE_INDEX, unit))
                end
            end
        end
    end
    fprintf('\n------------------------------\n')
end
