% Sept 27 2012: Created code to provide consistency and aesthetic nice-ness
% to the many histograms that Larry has been working on for the pursuit vs
% saccade signal project in FEF.  Goals are to make consistent and centered
% x-axis, 0.1 bin size, display mean on graph; stats/pval for mean, display
% sample size, std dev, other stuff.  JPM

% DEPENDENT FUNCTION
% 'prettyFigSGL' 


function polarhistogram_withSignif(column_of_values,bs,faceColor,lineColor,fig_name)

edges = deg2rad(0:bs:360);
% Plot histogram with white bars, range -1 to 
h = polarhistogram(deg2rad(column_of_values), 'BinEdges',edges,'FaceColor', faceColor,'FaceAlpha',0.5,'Normalization','pdf');
hold on

% Plot star at MEAN, 0.5 above max of ylim
polarplot(repmat(mean(h.Data),1,2), [0 1], '--', 'color', lineColor,'linewidth',2)
polarplot(repmat(median(h.Data),1,2), [0 1], ':', 'color', lineColor,'linewidth',2)

% title(fig_name,'Interpreter','tex') 
% xlabel(xlab) 
% ylabel(ylab)                

%prettyFigSGL
prettyFig

% Perform signrank test to determine if distribution is significantly
% different from 0 (Note the technical wrinkle that the signrank test used
% the median but we are displaying the mean; this is standard)
datamean = mean(column_of_values,'omitnan');
datamedian = median(column_of_values,'omitnan');
datastd = std(column_of_values,'omitnan');
[pvalue, hypothesis] = signrank(column_of_values);

legend(sprintf('n = %d', length(column_of_values)),sprintf('mean = %0.3f',datamean),sprintf('median = %0.3f', datamedian),'location','southoutside','orientation','horizontal')
title([fig_name, sprintf(' (p = %0.3f)',pvalue)],'fontsize',14)

% Plot values on figure
% if hypothesis == 1 % if null hypothesis is rejected/ statistically signif, use green font for pvalue
%     text (xlPos, max(ylims)*0.97, ['n = ', num2str(length(column_of_values))], 'fontsize', 15, 'color', 'k' ) % sample size
%     text (xlPos, max(ylims)*0.91, ['mean = ', sprintf('%0.3f', datamean)], 'fontsize', 15, 'color', 'k' ) % mean
%     text (xlPos, max(ylims)*0.85, ['std = ', sprintf('%0.3f', datastd)], 'fontsize', 15, 'color', 'k' ) % std
%     text (xlPos, max(ylims)*0.79, ['p = ', sprintf('%0.3f', pvalue)], 'fontsize', 15, 'color', 'k', 'fontweight', 'bold' ) % pvalue
% else % If not significant, pvalue in black
%     text (xlPos, max(ylims)*0.97, ['n = ', num2str(length(column_of_values))], 'fontsize', 15, 'color', 'k' ) % sample size
%     text (xlPos, max(ylims)*0.91, ['mean = ', sprintf('%0.3f', datamean)], 'fontsize', 15, 'color', 'k' ) % mean
%     text (xlPos, max(ylims)*0.85, ['std = ', sprintf('%0.3f', datastd)], 'fontsize', 15, 'color', 'k' ) % std
%     text (xlPos, max(ylims)*0.79, ['p = ', sprintf('%0.3f', pvalue)], 'fontsize', 15, 'color', 'k' ) % pvalue
% end