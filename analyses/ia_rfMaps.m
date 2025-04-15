function ia_rfMaps(data_path,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data_path', @ischar);
    addParameter(p, 'FIG_PATH', @ischar);
    addParameter(p, 'IMEC', 0, @isnumeric);
    
    parse(p, data_path, varargin{:});
    data_path = p.Results.data_path;
    FIG_PATH = p.Results.FIG_PATH;
    IMEC = p.Results.IMEC;

    fprintf('\n------------------------------\n')
    load(data_path,'S');
    [parent_path, filename, ~] = fileparts(data_path);

    if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end
    
    TBL = S.rfmp1.data;
    
    for unit=1:height(S.kilosort(IMEC+1).clusters)
        if ~exist(fullfile(FIG_PATH, sprintf('imec%d_unit%03d.png', IMEC, unit)), 'file')
            [frs,bin_edges,xvals,yvals] = format_tableToRFMap(TBL, 'IMEC', IMEC, 'UNITS', unit);
        
            f2a = figure; %('Visible','off');
            f2a.Position = [100 100 1800 900];
            tl = heatMap_rfOverTime(frs{1},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',pix2deg(xvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)), 'Y_VALS',pix2deg(yvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)));
            
            if (IMEC)==0
                title(tl,sprintf('%s --- LEFT --- unit %d',filename, unit),'fontsize',20,'interpreter','none')
            else
                title(tl,sprintf('%s --- RIGHT --- unit %d',filename, unit),'fontsize',20,'interpreter','none')
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
