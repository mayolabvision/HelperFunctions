function pieChart_trialOutcomes(trialOutcomes, taskname, all_results)

    % Sort 'all_results' alphabetically for consistent color mapping
    sorted_all_results = categorical(sort(categories(all_results)));

    % Sort trialOutcomes based on the sorted_all_results order
    sorted_trialOutcomes = [];
    for s = 1:length(sorted_all_results)
        if ismember(sorted_all_results(s), trialOutcomes)
            sorted_trialOutcomes = [sorted_trialOutcomes; trialOutcomes(ismember(trialOutcomes, sorted_all_results(s)))];
        end
    end

    % Assign consistent colors to the categoricals in 'all_results'
    color_map = containers.Map; % Map for consistent color assignments
    for i = 1:length(sorted_all_results)
        category = sorted_all_results(i);
        switch char(category)
            case 'CORRECT'
                color_map(char(category)) = [0, 0.5, 0]; % Green
            case 'IGNORED'
                color_map(char(category)) = [0.5, 0.5, 0.5]; % Gray
            otherwise
                % Assign a shade of pink/red for other categories
                rng(sum(double(char(category)))); % Seed the random generator for consistent colors
                color_map(char(category)) = [0.8 + rand() * 0.2, rand() * 0.3, rand() * 0.3];
        end
    end

    % Replace underscores in 'trialOutcomes' for cleaner labels
    escapedLabels = strrep(cellstr(sorted_trialOutcomes), '_', '\_');

    % Create pie chart
    h = pie(categorical(escapedLabels));

    % Assign colors to the pie chart slices manually using a map for consistency
    for i = 1:length(h)
        if strcmp(h(i).Type, 'patch')
            % Directly assign the color from the map
            outcome = h(i).DisplayName;
            if isKey(color_map, outcome)
                h(i).FaceColor = color_map(outcome); % Apply corresponding color
            end
        end
    end

    % Add title
    title({'';'';sprintf('%s (%d correct / %d total trials)', taskname, sum(trialOutcomes=="CORRECT"), length(trialOutcomes))},'fontsize',14);
end