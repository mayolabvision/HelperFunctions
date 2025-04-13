function ia_trialOutcomes(data_path)
    load(data_path,'S');
    [parent_path, filename, ~] = fileparts(data_path);

    fig_path = fullfile(parent_path, 'figs', 'session_trialOutcomes.png');

    fn = fieldnames(S);
    tasks = fn(~ismember(fn, {'kilosort', 'recording_info'}));

    results = cell(length(tasks),1);
    for task=1:length(tasks)
        results{task} = S.(tasks{task}).data.result;
    end
    
    f1a = figure('Visible','off');
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
    print(f1a, fig_path, '-dpng', '-r300');
    fprintf(sprintf('\n----Trial outcome summary plot for %s COMPLETE----\n',filename))

    fprintf('\n------------------------------\n')
    
end