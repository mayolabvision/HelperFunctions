% Sept 27 2012: Created code to provide consistency and aesthetic nice-ness
% to the many histograms that Larry has been working on for the pursuit vs
% saccade signal project in FEF.  Goals are to make consistent and centered
% x-axis, 0.1 bin size, display mean on graph; stats/pval for mean, display
% sample size, std dev, other stuff.  JPM

% DEPENDENT FUNCTION
% 'prettyFigSGL' 


function histogram_withSignif(column_of_values,xlims,ylims,bs,xlPos,lineColor,faceColor,includeP)

% Plot histogram with white bars, range -1 to 
h = histogram(column_of_values, 'binwidth',bs,'linewidth', 2, 'facecolor', faceColor); % -1:binsize:1,
xlim(xlims)
ylim(ylims)
hold on
axis square

% Plot star at MEAN, 0.5 above max of ylim
plot (mean(column_of_values,'omitnan'), max(h.Values)+max(h.Values)/20, '*', 'markersize', 15, 'color', 'k')

% draw thin vertical dotted line at 0 for reference
%line([0,0], ylims, 'color', 'k', 'linewidth', 2, 'linestyle', ':')

% draw thick vertical solid line at mean (redundant with star above)
line([mean(column_of_values,'omitnan'), mean(column_of_values,'omitnan') ], ylims, 'color', lineColor, 'linewidth', 3 )

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

% Plot values on figure
if hypothesis == 1 % if null hypothesis is rejected/ statistically signif, use green font for pvalue
    text (xlPos, max(ylims)*0.97, ['n = ', num2str(length(column_of_values))], 'fontsize', 15, 'color', 'k' ) % sample size
    text (xlPos, max(ylims)*0.91, ['mean = ', sprintf('%0.3f', datamean)], 'fontsize', 15, 'color', 'k' ) % mean
    text (xlPos, max(ylims)*0.85, ['std = ', sprintf('%0.3f', datastd)], 'fontsize', 15, 'color', 'k' ) % std
    if includeP==1
        text (xlPos, max(ylims)*0.79, ['p = ', sprintf('%0.3f', pvalue)], 'fontsize', 15, 'color', 'k', 'fontweight', 'bold' ) % pvalue
    end
else % If not significant, pvalue in black
    text (xlPos, max(ylims)*0.97, ['n = ', num2str(length(column_of_values))], 'fontsize', 15, 'color', 'k' ) % sample size
    text (xlPos, max(ylims)*0.91, ['mean = ', sprintf('%0.3f', datamean)], 'fontsize', 15, 'color', 'k' ) % mean
    text (xlPos, max(ylims)*0.85, ['std = ', sprintf('%0.3f', datastd)], 'fontsize', 15, 'color', 'k' ) % std
    if includeP==1
        text (xlPos, max(ylims)*0.79, ['p = ', sprintf('%0.3f', pvalue)], 'fontsize', 15, 'color', 'k' ) % pvalue
    end
end