clear
clc

% Directory containing the .mat files
folderPath = '/Users/kendranoneman/Data/NP_DATA/param_sweep';  % <-- Change this to your folder
fileList = dir(fullfile(folderPath, '*.mat'));

allClustersTable = table();
for i = 1:length(fileList)
    filename = fileList(i).name;
    filepath = fullfile(folderPath, filename);
    fprintf('\n~~~%s~~~\n',filename)

    % Extract whitening_range from filename using regular expression
    tokens = regexp(filename, 'wh(\d+)', 'tokens');
    whitening_range = str2double(tokens{1}{1});

    % Extract whitening_range from filename using regular expression
    tokens = regexp(filename, 'ct(\d+)', 'tokens');
    ccg_thresh = str2double(tokens{1}{1});

    % Load the S struct from the file
    fileData = load(filepath, 'S');

    clusterTable = fileData.S.kilosort.clusters;
    clusterTable.whitening_range = repmat(whitening_range, height(clusterTable), 1);
    clusterTable.ccg_thresh = repmat(ccg_thresh, height(clusterTable), 1);
    clusterTable.num_clusts = repmat(height(clusterTable), height(clusterTable),1);
    clusterTable.perc_good = repmat((sum(clusterTable.KSLabel_clusters=='good')/height(clusterTable))*100,height(clusterTable),1);


    allClustersTable = [allClustersTable; clusterTable];
end

save(fullfile(folderPath,'allClustersTable.mat'), 'allClustersTable');

%% PLOT SWEEP DATA
f1 = figure;
f1.Position = [100 100 900 1000];

% Define metrics to plot
metrics = {'num_clusts','perc_good', 'snr', 'fr_hz', 'vis_sel_dir', 'sac_sel_dir'};
metric_names = {'# clusters', '% "good" clusters', 'SNR', 'FR [Hz]', 'dir sel (vis)', 'dir sel (mot)'};
numMetrics = numel(metrics);

% Use tiledlayout for shared x-axis
tl = tiledlayout(numMetrics, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xlabelShared = ["Whitening Range (# of nearby channels)", "CCG Threshold (% RF violations allowed)"];

allAxes = gobjects(numMetrics, 2);  % to store axes handles for linking

% Group by whitening_range or ccg_thresh
for m = 1:2
    if m==1
        tbl = allClustersTable(allClustersTable.ccg_thresh==25,:);
        [G, groupIDs] = findgroups(tbl.whitening_range);
    else
        tbl = allClustersTable(allClustersTable.whitening_range==14,:);
        [G, groupIDs] = findgroups(tbl.ccg_thresh);
    end
    x = groupIDs;

    ax_cnt = 1;
    for i = 1:numMetrics
        ax_cnt = ax_cnt + 1;
        metric = metrics{i};
        metric_name = metric_names{i};

        % Compute mean and SEM
        meanVals = splitapply(@nanmean, tbl.(metric), G);
        semVals  = splitapply(@(x) nanstd(x) / sqrt(sum(~isnan(x))), tbl.(metric), G);

        if m==1
            ax1(ax_cnt) = nexttile((i - 1) * 2 + m);  % column-major order
        else
            ax2(ax_cnt) = nexttile((i - 1) * 2 + m);  % column-major order
        end
        errorbar(x, meanVals, semVals, 'ko-', 'LineWidth', 1.5);

        if i == numMetrics
            xlabel(xlabelShared(m));
        end
        if m == 1
            ylabel(metric_name, 'Interpreter', 'none');
            xlim([1,33])
            xticks([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32])
            xticklabels({"2","","6","","10","","14","","18","","22","","26","","30",""})
        else
            xlim([0,26])
            xticks([1,5,9,13,17,21,25])
            xticklabels({"1","5","9","13","17","21","25"})
        end
        prettyFig;
    end
end

linkaxes(ax1,'x');
linkaxes(ax2,'x');
for row = 1:numMetrics
    linkaxes([ax1(row+1),ax2(row+1)],'y')
end

print(f1, fullfile(folderPath, 'SWEEP_RESULTS.png'), '-dpng', '-r200');
