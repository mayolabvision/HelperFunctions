function pieChart_trialOutcomes(trialOutcomes)

escapedLabels = strrep(cellstr(trialOutcomes), '_', '\_');

piechart(categorical(escapedLabels), FaceAlpha=0.85); % Use the modified labels

newcolors = [65 65 65
             75 75 75
             0 0 0]./255;
         
colororder(newcolors)

title(sprintf('Trial Outcomes (N = %d total trials)', length(trialOutcomes)));

end