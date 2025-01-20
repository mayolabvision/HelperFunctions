function pieChart_trialOutcomes(trialOutcomes)

escapedLabels = strrep(cellstr(trialOutcomes), '_', '\_');

piechart(categorical(escapedLabels)); % Use the modified labels

title(sprintf('Trial Outcomes (N = %d total trials)', length(trialOutcomes)));

end