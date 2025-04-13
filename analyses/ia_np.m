%% 
OUT_PATH='/Volumes/home/DATA';

EXPERIMENTER = 'kendra';
MONKEY = 'scrappy';
SESSION = '0129a';

filename = sprintf('%s_%s_%s_g0',EXPERIMENTER,MONKEY,SESSION);
data_path = fullfile(OUT_PATH,filename,[filename,'.mat']);

% load(sprintf('/Volumes/home/DATA/%s/%s.mat',filename,filename))

% lr_colors = {[123, 44, 191]./255; [42, 157, 143]./255};
% lr_med = {[189,149,223]./255; [127,196,187]./255};
% lr_dark = {[86,30,133]./255; [33,125,114]./255};
% lr_shades = {[228,212,242]./255; [212,235,232]./255};
% 
% lr_comp_colors = {[191,44,186]./255; [42,114,157]./255};
% lr_comp_shades = {[242,212,241]./255; [212,226,235]./255};
% lr_comp_med = {[216,128,213]./255; [127,170,196]./255};
% lr_comp_dark = {[152,35,148]./255; [33,91,125]./255};

%% ------------------------------------------- 1. Overall Performance & BEHAVIOR plots -------------------------------------------
filename = sprintf('%s_%s_%s_g0',EXPERIMENTER,MONKEY,SESSION);
TABLE_PATH = fullfile(OUT_PATH,filename,[filename,'.mat']);

ia_trialOutcomes(TABLE_PATH)

%% ------------------------------------------- 2. RF MAPPING -------------------------------------------
FILENAME = 'kendra_scrappy_0129a_g0';
DATA_PATH = fullfile(OUT_PATH,FILENAME,[FILENAME,'.mat']);
IMEC = 0;

ia_rfMaps(DATA_PATH,'IMEC',IMEC)

%% ------------------------------------------- 2. MG SACCADES -------------------------------------------
FILENAME = 'kendra_scrappy_0129a_g0';
DATA_PATH = fullfile(OUT_PATH,FILENAME,[FILENAME,'.mat']);
IMEC = 0;
ALIGN = 'stim';

ia_mdirRasters(DATA_PATH,'IMEC',IMEC,'ALIGN',ALIGN)

%% ------------------------------------------- 2. MG SACCADES -------------------------------------------
FILENAME = 'kendra_scrappy_0129a_g0';
DATA_PATH = fullfile(OUT_PATH,FILENAME,[FILENAME,'.mat']);
IMEC = 0;
ALIGN = 'sacc';

ia_mdirRasters(DATA_PATH,'IMEC',IMEC,'ALIGN',ALIGN)

%% fdsf    

for unit = 1:length(good_chans)
    chan_name   =  S.channels.mapped_name{good_chans(unit)};
    chan_depth  =  S.channels.depth_mm(good_chans(unit));

    for speed = 1:length(pursuitSpeeds)
        f3f = figure('Visible','off');
        f3f.Position = [100 100 1800 900];
        
        angles = sort(unique(T.angle))';
        angle_order = [6,3,2,1,4,7,8,9];
        y_lims = []; % Store y-axis limits
        frs_perAng = cell(length(angles),1);
        for ang = 1:length(angles)
            these_trls = T(T.angle==angles(ang) & T.pursuitSpeed==pursuitSpeeds(speed) & T.jump==jump & T.csFlag==csFlag,:);
        
            sptimes = cellfun(@(q,z,r) q(z)-r, these_trls.spiketimes(:,good_chans(unit)), cellfun(@(q) q>GAMMA, these_trls.net_labels(:,good_chans(unit)), 'uni', 0), num2cell(these_trls.PURSUIT_TARG), 'uni', 0);
        
            subplot(3,3,angle_order(ang))
        
            line_color = [0,0,0]./255; sem_shade = [200,200,200]./255;
            raster_sdf(sptimes', [50, 500], 10, 'line_color', line_color, 'sem_shade', sem_shade)
        
            yyaxis left;
            ax = gca;
            y_lims = [y_lims; ax.YLim];
        
            frs_perAng{ang} = cellfun(@(q) (sum(q>=0 & q <500)*(1000/500)), sptimes, 'uni', 1);
        
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
            rhoPst = [rhoPst; nanmean(permutedStimrate)];
        end
        
        sorted_rhoPst=sort(rhoPst);
        rhoLst = sorted_rhoPst(shuffles*.05,:); % 95% lower confidence interval
        rhoUst = sorted_rhoPst(shuffles-(shuffles*.05),:); % 95% upper confidence interval
        
        % calculate tuning preferences
        %theta = 0:360/length(a.CND):360; theta(end)=[];
        theta = 0:45:315;
        [visds, visdp] = tuningbias(theta,nanmean(stimrate));
        
        subplot(3,3,5)
        rho = nanmean(stimrate);
        dst = sprintf('%0.2f',visds);
        dpt = sprintf('%0.2f',visdp);
        polarplot(deg2rad([theta 0]),[rho rho(1)],'ko-',...
            'markerfacecolor','k','linewidth',3)
        hold on
        polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'go--','LineWidth',2);
        polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'go--','LineWidth',2);   
        
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
        
        han=axes(f3f,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        xlabel(han,{'';'time aligned to target motion onset (ms)'},'fontsize',16);
        title(han,{sprintf('%s_purs',filename); sprintf('%s (ripChan = %d, depth = %2.3f mm)',chan_name, S.channels.ripChan_num(good_chans(unit)), chan_depth); sprintf('Pure Pursuit Trials (Speed = %d deg/s, Jump = %d)', pursuitSpeeds(speed), jump)},'fontsize',18,'interpreter','none')
    
        print(f3f, fullfile(fig_path, sprintf('%s-purePursuit-s%.2d.png', chan_name, pursuitSpeeds(speed))), '-dpng', '-r200');
        fprintf(sprintf('\n----Unit %.2d complete----',good_chans(unit)))
    
    end
end
fprintf('\n----------------------\n')









