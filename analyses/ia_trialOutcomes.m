function ia_trialOutcomes(data_path,varargin)
    p = inputParser;
    addRequired(p, 'data_path', @ischar);
    addParameter(p, 'FIG_PATH', [], @ischar);
    
    parse(p, data_path, varargin{:});
    data_path = p.Results.data_path;
    FIG_PATH = p.Results.FIG_PATH;

    fprintf('\n------------------------------\n')
    load(data_path,'S');
    [parent_path, filename, ~] = fileparts(data_path);
    fprintf(sprintf('\n----Data loaded for %s----\n',filename))
   
    if isempty(FIG_PATH)
        FIG_PATH = fullfile(parent_path, 'figs');
    end
    if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end   

    fn = fieldnames(S);
    tasks = fn(~ismember(fn, {'kilosort', 'recording_info'}));

    results = cell(length(tasks),1);
    for task=1:length(tasks)
        results{task} = S.(tasks{task}).tbl.result;
    end
    
    f1a = figure; %('Visible','off');
    f1a.Position = [100 100 500*length(tasks) 400];
    tl = tiledlayout(1,length(tasks));
    tl.TileSpacing = 'tight';
    tl.Padding = 'compact';
    
    set(gcf,'color','w')
    for ii = 1:length(tasks)
        nexttile
        pieChart_trialOutcomes(results{ii}, tasks{ii}, vertcat(results{:}));
    end
    title(tl,sprintf('%s',filename),'fontsize',18,'interpreter','none')
    print(f1a, fullfile(FIG_PATH,'session_trialOutcomes.png'), '-dpng', '-r300');
    fprintf(sprintf('\n----Trial outcome summary plot for %s COMPLETE----\n',filename))

    fprintf('\n------------------------------\n')
    
end
