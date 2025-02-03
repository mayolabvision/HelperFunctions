function tl = polarPlot_pursDistByCondition(T,propType)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    tl = tiledlayout(1,2);
    tl.TileSpacing = 'loose';
    tl.Padding = 'tight';

    % percent correct by condition
    angs = sort(unique(T.angle)); angs = [angs; angs(1)];
    colors = cellfun(@(q) radial_colormap(q), num2cell(angs), 'uni', 0);
    
    jumps = sort(unique(T.jump));
    speeds = sort(unique(T.pursuitSpeed));
    perDists = zeros(length(angs),length(jumps),length(speeds));
    for a=1:length(angs)
        for j=1:length(jumps)
            for s=1:length(speeds)
                if isequal(propType,'perCorrect')
                    perDists(a,j,s) = 100 * sum(T.result=='CORRECT' & T.angle==angs(a) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s))/sum(T.angle==angs(a) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s));
                    rlim = [0 100];
                    rticks = [0 20 40 60 80 100];
                    rlabelPos = 40;
                    rLabel = '% correct';
                elseif isequal(propType,'perPure')
                    perDists(a,j,s) = 100* sum(T.csFlag==0 & T.result=='CORRECT' & T.angle==angs(a) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s)) / sum(T.result=='CORRECT' & T.angle==angs(a) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s));
                    rlim = [0 100];
                    rticks = [0 20 40 60 80 100];
                    rlabelPos = 40;
                    rLabel = '% pure';
                elseif isequal(propType,'pursLatency')
                    perDists(a,j,s) = mean(T.rxnTime(T.result=='CORRECT' & T.msFlag==0 & T.angle==angs(a) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s)));
                    rlim = [50 175];
                    rticks = [50 75 100 125 150 175 200 225 250];
                    rlabelPos = 150;
                    rLabel = 'latency';
                end
            end
        end
    end
    
    angs_rad = angs * (pi / 180);
    ls = {'-','--',':'};
    
    for s=1:length(speeds)
        nexttile
        l = gobjects(1, length(jumps)); % Preallocate array of graphics objects
        for j=1:length(jumps)
            this_line = perDists(:,j,s);
    
            l(j) = polarplot(angs_rad,this_line, 'Color', 'black', 'LineStyle', ls{j}, 'LineWidth', 1.5);
            hold on;
            % Plot each point with its specific color
            for a = 1:length(angs)
                polarplot(angs_rad(a), this_line(a), 'o', ...
                    'MarkerSize', 10, 'MarkerFaceColor', colors{a}, 'MarkerEdgeColor', colors{a});
            end
        end
        ax = gca; % Get the polar axes
        ax.ThetaTick = 0:30:180; % Set angular ticks (degrees)
        ax.RLim = rlim; % Set radial limits
        ax.RTick = rticks; % Set radial ticks
    
        if s==1
            ll = legend(l, string(jumps)');  % Create the legend
            title(ll, 'Jump Values');
        end
    
        % Add labels for the axes
        % Radial label

        % text(0, rlabelPos, rLabel, 'HorizontalAlignment', 'center', ...
        %     'VerticalAlignment', 'bottom', 'FontSize', 12);
    
        title(sprintf('speed = %d deg/s',speeds(s)),'fontsize',16)
        prettyFig;
    end
end