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

    if IMEC==0
        line_color = [123, 44, 191]./255;
        sem_shade = [228,212,242]./255;
    else
        line_color = [42, 157, 143]./255;
        sem_shade = [212,235,232]./255;
    end
    
    if isfield(S, 'mdir2')
        T = [S.mdir1.data; S.mdir2.data];
    else
        T = S.mdir1.data;
    end
    T = T(T.result=='CORRECT',:);
    
    angles = sort(unique(T.angle))';
    angle_order = [6,3,2,1,4,7,8,9];

    if isnan(JOB_ID)
        units = 1:height(S.kilosort(IMEC+1).clusters);
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
        if ~exist(fullfile(fig_path, sprintf('imec%d_unit%03d.png', IMEC, unit)), 'file')
            f3a = figure; %('Visible','off');
            f3a.Position = [100 100 1800 900];
        
            y_lims = []; % Store y-axis limits
            frs_perAng = cell(length(angles),1);
            for ang = 1:length(angles)
                these_trls = T(T.angle==angles(ang),:);
                
                if isequal(ALIGN,'stim')
                    sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(imec_name), 'uni', 0), these_trls.TARG_ON, 'uni', 0);
                elseif isequal(ALIGN,'sacc')
                    sptimes = cellfun(@(w,v) w-v(1), cellfun(@(q) q{(unit+1)}, these_trls.(imec_name), 'uni', 0), num2cell(these_trls.SACCADE), 'uni', 0);
                end

                subplot(3,3,angle_order(ang))
                raster_sdf(sptimes', 'TIME_WINDOW', X_LIMITS, 'LINE_COLOR', line_color, 'SEM_SHADE', sem_shade)
        
                yyaxis left;
                ax = gca;
                y_lims = [y_lims; ax.YLim];
            
                frs_perAng{ang} = cellfun(@(q) (sum(q>=FR_WIN(1) & q <FR_WIN(2))*(1000/500)), sptimes, 'uni', 1);
        
            end
        
            % Pad each cell with NaNs to match maxLength
            maxLength = max(cellfun(@numel, frs_perAng));
            frs_perAng = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_perAng, 'UniformOutput', false);
            
            stimrate = vertcat(frs_perAng{:})';
            
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
            
            subplot(3,3,5)
            rho = mean(stimrate,'omitnan');
            dst = sprintf('%0.2f',visds);
            dpt = sprintf('%0.2f',visdp);
            polarplot(deg2rad([theta 0]),[rho rho(1)],'ko-',...
                'markerfacecolor','k','linewidth',3)
            hold on
            polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'o--','LineWidth',2,'Color',line_color);
            polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'o--','LineWidth',2,'Color',line_color);   
            
            h1=polarplot(deg2rad(visdp),max(rho),'k^'); % Plot black triangle at best dir
            txd1 = get(h1,'ThetaData');
            tyd1 = get(h1,'RData');
            set(h1,'ThetaData',txd1(1),'RData',tyd1(1));
            set(h1,'markersize',10,'markerfacecolor','k'); hold off;
            title(['VisDir: ',dpt,', Sel: ',dst]);
            prettyFig
            
            
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
        
            print(f3a, fullfile(fig_path, sprintf('imec%d_unit%03d.png', (IMEC), unit)), '-dpng', '-r200');
            fprintf(sprintf('\n----IMEC %d, Unit %.3d COMPLETE----',IMEC, unit))
        else
            fprintf(sprintf('\n----IMEC %d, Unit %.3d exists----',IMEC, unit))
        end

    end
    fprintf('\n------------------------------\n')
end
