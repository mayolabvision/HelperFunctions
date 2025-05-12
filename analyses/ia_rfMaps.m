function ia_rfMaps(data_path,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data_path', @ischar);
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'IMEC', 0, @isnumeric);
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    addParameter(p, 'N_CHUNKS', NaN, @isnumeric);
    
    parse(p, data_path, varargin{:});
    data_path = p.Results.data_path;
    FIG_PATH = p.Results.FIG_PATH;
    IMEC = p.Results.IMEC;
    JOB_ID = p.Results.JOB_ID;
    N_CHUNKS = p.Results.N_CHUNKS;

    fprintf('\n------------------------------\n')
    load(data_path,'S');
    [parent_path, filename, ~] = fileparts(data_path);
    fprintf(sprintf('\n----Data loaded for %s----\n',filename))

    if isempty(FIG_PATH)
        FIG_PATH = fullfile(parent_path, 'figs', 'rfmp', 'unit_heatmaps');
    end
    if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end    

    units = S.kilosort(IMEC+1).clusters.cluster_id;
    chans = S.kilosort(IMEC+1).clusters.channel_id;
    snrs  = S.kilosort(IMEC+1).clusters.snr;
    contams = S.kilosort(IMEC+1).clusters.ContamPct;
    kslabs = S.kilosort(IMEC+1).clusters.KSLabel_cc;

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
        snr = snrs(ids);
        contams = contams(ids);
        kslabs = kslabs(ids);
        
    end
    
    % Find rfmp or rfMapping fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'rfmp', 'rfMapping'}, 'IgnoreCase', true));

    T = []; 
    for mm = 1:numel(matchingFields)
        T = [T; S.(matchingFields{mm}).tbl];
    end

    for u=1:length(units)
        unit = units(u); 
        if ~exist(fullfile(FIG_PATH, sprintf('imec%d_unit%04d_chan%03d.png', IMEC, unit, chans(u))), 'file')
            [frs,bin_edges,xvals,yvals] = format_tableToRFMap(T, 'IMEC', IMEC, 'UNITS', (unit+1));
        
            f2a = figure('Visible','off');
            f2a.Position = [100 100 1800 900];
            tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',xvals, 'Y_VALS',yvals);
            
            if (IMEC)==0
                title(tl,sprintf('%s --- LEFT --- cluster %d (channel %d)',filename, unit, chans(u)),'fontsize',20,'interpreter','none')
            else
                title(tl,sprintf('%s --- RIGHT --- cluster %d (channel %d)',filename, unit, chans(u)),'fontsize',20,'interpreter','none')
            end
            subtitle(tl, sprintf('ks_label = %s, snr = %.2f uV, contam_pct = %.1f%%', kslabs{u}, snrs(u), contams(u)),'fontsize',16,'interpreter','none')
            
            annotation('textbox', [0.75 0.89 0.2 0.1], ... % [x y w h] in normalized figure units
                'String', sprintf('N = %d repeats', min(min(min(cellfun(@length, frs{1}))))), ...
                'FontSize', 16, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'right');
            

            print(f2a, fullfile(FIG_PATH, sprintf('imec%d_unit%04d_chan%03d.png', IMEC, unit, chans(u))), '-dpng', '-r200');
            fprintf(sprintf('\n----IMEC %d, Unit %.4d COMPLETE----',IMEC, unit))
        else
            fprintf(sprintf('\n----IMEC %d, Unit %.4d exists----',IMEC, unit))
        end
    end
    fprintf('\n------------------------------\n')
end
