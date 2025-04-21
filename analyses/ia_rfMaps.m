function ia_rfMaps(data_path,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data_path', @ischar);
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'IMEC', 0, @isnumeric);
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    
    parse(p, data_path, varargin{:});
    data_path = p.Results.data_path;
    FIG_PATH = p.Results.FIG_PATH;
    IMEC = p.Results.IMEC;
    JOB_ID = p.Results.JOB_ID;

    fprintf('\n------------------------------\n')
    load(data_path,'S');
    [parent_path, filename, ~] = fileparts(data_path);

    if isempty(FIG_PATH)
        FIG_PATH = fullfile(parent_path, 'figs', 'rfmp', 'unit_heatmaps');
    end
    if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end    

    if isnan(JOB_ID)
        units = 1:height(S.kilosort(IMEC+1).clusters);
    else
        all_units = 1:height(S.kilosort(IMEC+1).clusters);
        % Split into 50 chunks as a cell array
        chunks = arrayfun(@(i) all_units(...
            floor((i-1)*numel(all_units)/50)+1 : ...
            floor(i*numel(all_units)/50)), ...
            1:50, 'UniformOutput', false);
        units = (chunks{(JOB_ID+1)}) - 1;
    end
    
    for u=1:length(units)
        unit = units(u);
        if ~exist(fullfile(FIG_PATH, sprintf('imec%d_unit%03d.png', IMEC, unit)), 'file')
            [frs,bin_edges,xvals,yvals] = format_tableToRFMap(S.rfmp1.data, 'IMEC', IMEC, 'UNITS', (unit+1));
        
            f2a = figure('Visible','off');
            f2a.Position = [100 100 1800 900];
            tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',pix2deg(xvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)), 'Y_VALS',pix2deg(yvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)));
            
            if (IMEC)==0
                title(tl,sprintf('%s --- LEFT --- cluster %d',filename, unit),'fontsize',20,'interpreter','none')
            else
                title(tl,sprintf('%s --- RIGHT --- cluster %d',filename, unit),'fontsize',20,'interpreter','none')
            end
            %subtitle(tl,sprintf('%s (ripChan = %d, depth = %2.3f mm)',chan_name, S.channels.ripChan_num(good_chans(unit)), chan_depth),'fontsize',16,'interpreter','none')
            
            print(f2a, fullfile(FIG_PATH, sprintf('imec%d_unit%03d.png', IMEC, unit)), '-dpng', '-r200');
            fprintf(sprintf('\n----IMEC %d, Unit %.3d COMPLETE----',IMEC, unit))
        else
            fprintf(sprintf('\n----IMEC %d, Unit %.3d exists----',IMEC, unit))
        end
    end
    fprintf('\n------------------------------\n')
end
