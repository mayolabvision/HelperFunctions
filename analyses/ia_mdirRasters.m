function ia_mdirRasters(data_path,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    p = inputParser;
    addRequired(p, 'data_path', @ischar);
    addParameter(p, 'IMEC', 0, @isnumeric);
    addParameter(p, 'FIG_PATH', [], @ischar);
    addParameter(p, 'ALIGN', 'stim', @ischar);
    addParameter(p, 'X_LIMITS', [-300 500], @isnumeric)
    addParameter(p, 'JOB_ID', NaN, @isnumeric);
    
    parse(p, data_path, varargin{:});
    data_path = p.Results.data_path;
    IMEC = p.Results.IMEC;
    FIG_PATH = p.Results.FIG_PATH;
    ALIGN = p.Results.ALIGN;
    X_LIMITS = p.Results.X_LIMITS;
    JOB_ID = p.Results.JOB_ID;

    load(data_path,'S');
    [parent_path, filename, ~] = fileparts(data_path);

    if isempty(FIG_PATH)
        FIG_PATH = fullfile(parent_path, 'figs', 'mdir', 'unit_rasters');
    end

    if isequal(ALIGN,'stim')
        FR_WIN = [50,150];
        fig_path = fullfile(FIG_PATH, 'stim_aligned');
        xlab = 'time aligned to target onset (ms)';
    elseif isequal(ALIGN,'sacc')
        FR_WIN = [-50,50];
        fig_path = fullfile(FIG_PATH, 'sacc_aligned');
        xlab = 'time aligned to saccade onset (ms)';
    end
    if ~exist(fig_path, 'dir'), mkdir(fig_path); end    


    % Find rfmp or rfMapping fields
    fields = fieldnames(S);
    matchingFields = fields(contains(fields, {'mdir', 'dirmem'}, 'IgnoreCase', true));
    
    T = []; 
    for mm = 1:numel(matchingFields)
        T = [T; S.(matchingFields{mm}).data];
    end

    T = T(T.result=='CORRECT',:);
    
    angles = sort(unique(T.angle))';
    angle_order = [6,3,2,1,4,7,8,9];
    distances = sort(unique(T.distance));

    if IMEC==0 % purples/pinks
        line_color = {[123,44,191]./255; [230,34,172]./255};
        tick_color = {[98,35,152]./255; [184,27,137]./255};
        sem_shade = {[228,212,242]./255; [250,210,238]./255};
    else % greens/blues
        line_color = {[42,157,143]./255; [42,114,157]./255};
        tick_color = {[25,94,85]./255; [25,68,94]./255};
        sem_shade = {[212,235,232]./255; [212,226,235]./255};
    end
    

    if isnan(JOB_ID)
        units = (1:height(S.kilosort(IMEC+1).clusters)) - 1;
    else
        all_units = 1:height(S.kilosort(IMEC+1).clusters);
        % Split into 50 chunks as a cell array
        n_chunks = 100;
        chunks = arrayfun(@(i) all_units(...
            floor((i-1)*numel(all_units)/n_chunks)+1 : ...
            floor(i*numel(all_units)/n_chunks)), ...
            1:n_chunks, 'UniformOutput', false);
        units = (chunks{(JOB_ID+1)}) - 1;
    end
    
    imec_name = ['spiketimes_imec' num2str(IMEC)];
    for u=1:length(units)
        unit = units(u);
        if ~exist(fullfile(fig_path, sprintf('imec%d_unit%04d.png', IMEC, unit)), 'file')
            f3a = figure('Visible','off');
            f3a.Position = [100 100 1800 900];
        
            y_lims = []; % Store y-axis limits
            frs_perAng = cell(length(angles),length(distances));
            for ang = 1:length(angles)
                these_trls = T(T.angle==angles(ang),:);

                these_trls.trialName = categorical(regexprep( cellstr(these_trls.trialName), 'mdir1\.(\d+)', 'mdir1.${sprintf(''%04d'', str2double($1))}'));
                these_trls = sortrows(these_trls, {'distance', 'trialName'});

                if isequal(ALIGN,'stim')
                    sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(imec_name), 'uni', 0), these_trls.TARG_ON, 'uni', 0);
                elseif isequal(ALIGN,'sacc')
                    sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(imec_name), 'uni', 0), num2cell(these_trls.SACCADE), 'uni', 0);
                end

                subplot(3,3,angle_order(ang))

                if numel(distances)==1
                    raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color{1}, 'SEM_SHADE', sem_shade{1})
                    frs_perAng{ang} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/(FR_WIN(2)-FR_WIN(1)))), sptimes, 'uni', 1);
                else
                    [line_colors, tick_colors, sem_shades] = deal(cell(height(these_trls),1)); 
                    for dd = 1:numel(distances)
                        line_colors(these_trls.distance==distances(dd)) = line_color(dd);
                        tick_colors(these_trls.distance==distances(dd)) = tick_color(dd);
                        sem_shades(these_trls.distance==distances(dd)) = sem_shade(dd);

                        frs_perAng{ang,dd} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/(FR_WIN(2)-FR_WIN(1)))), sptimes(these_trls.distance==distances(dd)), 'uni', 1);
                    end

                    raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_colors, 'SEM_SHADE', sem_shades, 'TICK_COLOR', tick_colors, 'FR_WINDOW', FR_WIN)
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
                str_title{dd} = sprintf('%d deg -- VisDir: %0.2f, Sel: %0.2f', distances(dd), visdp, visds);
            end 

            title(str_title);
               
            legend_labels = arrayfun(@(d) sprintf('%d deg', distances(d)), 1:length(distances), 'UniformOutput', false);
            legend(h1, legend_labels, 'Location', 'best');
            prettyFig;

            % Find the global y-axis limits
            global_y_lim = [min(y_lims(:,1)), max(y_lims(:,2))];
            
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
        
            if (IMEC)==0
                title(han,sprintf('%s --- LEFT --- cluster %d',filename,unit),'fontsize',20,'interpreter','none')
            else
                title(han,sprintf('%s --- RIGHT --- cluster %d',filename,unit),'fontsize',20,'interpreter','none')
            end
        
            print(f3a, fullfile(fig_path, sprintf('imec%d_unit%04d.png', (IMEC), unit)), '-dpng', '-r200');
            fprintf(sprintf('\n----IMEC %d, Unit %.4d COMPLETE----',IMEC, unit))
        else
            fprintf(sprintf('\n----IMEC %d, Unit %.4d exists----',IMEC, unit))
        end
    end
    fprintf('\n------------------------------\n')
end
