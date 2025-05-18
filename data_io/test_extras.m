clear
clc

data_path = '/Users/kendranoneman/Data/NP_DATA/kendra_scrappy_0124a_g0_unleashed.mat';
fig_path = '/Users/kendranoneman/Data/NP_DATA/kendra_scrappy_0124a_g0/figs';

%data_path = '/Users/kendranoneman/Data/NP_DATA/Ya_250506_s390_g0_unleashed.mat';

load(data_path,'S');

%% VMI

fields = fieldnames(S);
matchingFields1 = fields(contains(fields, {'dirmem', 'mdir'}, 'IgnoreCase', true));

if ~isempty(matchingFields1)
    Tmdir = []; 
    for mm = 1:numel(matchingFields1)
        Tmdir = [Tmdir; S.(matchingFields1{mm}).tbl];
    end
    
    Tmdir = Tmdir(Tmdir.result=='CORRECT',:);

    for imec = 1:size(S.kilosort,1)
        % VISUAL
        [vis_sel_dir, vis_pref_dir, ~, ~, frs_perAng_vis] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[50,150],'ALIGN_TO','stim','IMEC',imec-1);
        S.kilosort(imec).clusters.vis_sel_dir = vis_sel_dir;
        S.kilosort(imec).clusters.vis_pref_dir = vis_pref_dir;

        % MOTOR
        [sac_sel_dir, sac_pref_dir, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','IMEC',imec-1);
        S.kilosort(imec).clusters.sac_sel_dir = sac_sel_dir;
        S.kilosort(imec).clusters.sac_pref_dir = sac_pref_dir;

        % VMI
        VMI_per_unit = zeros(size(frs_perAng_sac,2),1);
        for unit = 1:size(frs_perAng_sac,2)
            visFR = frs_perAng_vis(:,unit);
            sacFR = frs_perAng_sac(:,unit);
        
            VMI_per_unit(unit) = (mean(vertcat(visFR{:})) - mean(vertcat(sacFR{:})))/(mean(vertcat(visFR{:})) + mean(vertcat(sacFR{:})));
        end
        S.kilosort(imec).clusters.VMI = VMI_per_unit;
    end
end

%% SPI

fields = fieldnames(S);
matchingFields2 = fields(contains(fields, {'pursuit', 'purs'}, 'IgnoreCase', true));

if ~isempty(matchingFields2)
    Tpurs = []; 
    for mm = 1:numel(matchingFields2)
        Tpurs = [Tpurs; S.(matchingFields2{mm}).tbl];
    end
    Tpurs = Tpurs(Tpurs.result=='CORRECT' & Tpurs.jump==-1 & Tpurs.pursType=='pure' & (isnan(Tpurs.msOffset) | Tpurs.msOffset<0),:);

    for imec = 1:size(S.kilosort,1)
        % MOTOR (PURSUIT)
        [pur_sel_dir, pur_pref_dir, ~, ~, frs_perAng_pur] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','IMEC',imec-1);
        S.kilosort(imec).clusters.pur_sel_dir = pur_sel_dir;
        S.kilosort(imec).clusters.pur_pref_dir = pur_pref_dir;

    end
end

%%
if ~isempty(matchingFields1) & ~isempty(matchingFields2)
    for imec = 1:size(S.kilosort,1)
        % MOTOR (SACCADE)
        [~, ~, ~, ~, frs_perAng_sac] = calculate_direction_tuning_from_tbl(Tmdir,'FR_WIN',[-50,50],'ALIGN_TO','sacc','IMEC',imec-1);

        % MOTOR (PURSUIT)
        [~, ~, ~, ~, frs_perAng_pur] = calculate_direction_tuning_from_tbl(Tpurs,'FR_WIN',[-50,50],'ALIGN_TO','purs','IMEC',imec-1);

        % VMI
        SPI_per_unit = zeros(size(frs_perAng_sac,2),1);
        for unit = 1:size(frs_perAng_sac,2)
            sacFR = frs_perAng_sac(:,unit);
            purFR = frs_perAng_pur(:,unit);
        
            SPI_per_unit(unit) = (mean(vertcat(sacFR{:})) - mean(vertcat(purFR{:})))/(mean(vertcat(sacFR{:})) + mean(vertcat(purFR{:})));
        end
        S.kilosort(imec).clusters.SPI = SPI_per_unit;
    end
end

%% Behavior

ia_trialOutcomes(data_path,'FIG_PATH',fig_path);

%% PURSUIT BEHAVIOR

% Pure pursuit only
f1 = figure; %('Visible','off');
f1.Position = [100 100 1500 450*numel(unique(tbl.pursuitSpeed))];

[tl,num_pure] = eyeTraces_pursSplitByConditions(tbl,unique(tbl.CROSSING_TIME(~isnan(tbl.CROSSING_TIME))),'PREINT',PURS_PREINT,'POSTINT',PURS_POSTINT,'MS_THRESH',MS_THRESH,'PURE_ONLY',true);

title(tl,sprintf('%s_purs',session),'fontsize',20,'interpreter','none')
subtitle(tl, sprintf('(# of total pure pursuit trials (w/ no MS) = %d, step-ramp duration = %d ms)',num_pure,unique(tbl.CROSSING_TIME(~isnan(tbl.CROSSING_TIME)))))

if ~isempty(FIG_PATH)
    print(f1, fullfile(FIG_PATH, 'eyeTraces_pureTrials.png'), '-dpng', '-r300');
end


%% SPI

fields = fieldnames(S);
matchingFields1 = fields(contains(fields, {'purs', 'pursuit'}, 'IgnoreCase', true));

Tpurs = []; 
for mm = 1:numel(matchingFields1)
    Tpurs = [Tpurs; S.(matchingFields1{mm}).tbl];
end

Tpurs = Tpurs(Tpurs.result=='CORRECT',:);
theta = sort(unique(Tpurs.angle))';

% PURSUIT
frs_perAng_targ = cell(length(theta),height(S.kilosort.clusters));
for a = 1:length(theta)
    this_ang = Tpurs(Tpurs.angle==theta(a),:);
    visFR = cellfun(@(w,v) cellfun(@(q) sum(q>=(v+50) & q<(v+150))/0.1, w, 'uni', 0), this_ang.spiketimes_imec0, num2cell(this_ang.PURSUIT_TARG_ON), 'uni', 0);
    frs_perAng_targ(a,:) = num2cell(cell2mat(vertcat(visFR{:})),1);
end

maxLength = max(cellfun(@numel, frs_perAng_targ)); maxLength = max(maxLength);
[purs_sel_dir, purs_pref_dir] = deal(zeros(size(frs_perAng_targ,2),1));
for unit = 1:size(frs_perAng_targ,2)
    frs_perAng2 = cellfun(@(x) [x; nan(maxLength - numel(x), 1)]', frs_perAng_targ(:,unit), 'UniformOutput', false);
                    
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
    [visds, visdp] = tuningbias(theta,mean(stimrate,'omitnan'));
    purs_sel_dir(unit) = visds; 
    purs_pref_dir(unit) = visdp;
end

% SPI
SPI_per_unit = zeros(size(frs_perAng_sac,2),1);
for unit = 1:size(frs_perAng_sac,2)
    sacFR = frs_perAng_sac(:,unit);
    purFR = frs_perAng_targ(:,unit);

    SPI_per_unit(unit) = (mean(vertcat(sacFR{:})) - mean(vertcat(purFR{:})))/(mean(vertcat(sacFR{:})) + mean(vertcat(purFR{:})));
end

%% RFMP

% Find rfmp or rfMapping fields
fields = fieldnames(S);
matchingFields1 = fields(contains(fields, {'rfmp', 'rfMapping'}, 'IgnoreCase', true));

Trfmp = []; 
for mm = 1:numel(matchingFields1)
    Trfmp = [Trfmp; S.(matchingFields1{mm}).tbl];
end

%%
% Parameters
dotRad = 45;                         % Radius of each circle
spacing = 2 * dotRad;               % Distance between centers of touching circles
xpos = -495:spacing:495;            % X positions
ypos = -495:spacing:495;            % Y positions
theta = linspace(0, 2*pi, 100);     % Angle for drawing circles

% Create figure
figure;
hold on;
axis equal;
axis([-540 540 -540 540]);         % Set axis limits a bit beyond 1080x1080 for padding
set(gca, 'Color', 'k');             % Optional: black background

% Draw circles at each (x,y) position
for x = xpos
    for y = ypos
        xc = x + dotRad * cos(theta);
        yc = y + dotRad * sin(theta);
        fill(xc, yc, 'w', 'EdgeColor', 'none');  % White circles with no border
    end
end

title('RFMapping task: 12x12 grid');
xlabel('x-position [pix]')
ylabel('y-position [pix]')
prettyFig;

%%
% Parameters for 24x24 grid
dotRad = 27;                        % Radius of each circle
spacing = 2 * dotRad;              % 44 pixels between centers
xpos = -513:spacing:513;           % 24 positions across 1080 px
ypos = -513:spacing:513;           % 24 positions down 1080 px
theta = linspace(0, 2*pi, 100);    % Circle angle resolution

% Create figure
figure;
hold on;
axis equal;
axis([-540 540 -540 540]);         % Slight padding for visibility
set(gca, 'Color', 'k');             % Optional: black background

% Draw the 24x24 grid of circles
for x = xpos
    for y = ypos
        xc = x + dotRad * cos(theta);
        yc = y + dotRad * sin(theta);
        fill(xc, yc, 'w', 'EdgeColor', 'none');  % White filled circles
    end
end

title('RFMapping task: 20x20 grid');
xlabel('x-position [pix]')
ylabel('y-position [pix]')
prettyFig;

