function S_new = unify_taskTables(S, tasknames)
    % Create a new struct to store the updated version of S
    S_new = S;

    % Loop through each task name in the tasknames cell array
    for t = 1:length(tasknames)
        task = tasknames{t};
        
        % Find fieldnames that contain the current task name
        fields = fieldnames(S);
        taskFields = fields(contains(fields, task));
        
        % Extract all column names from the tables for the current task
        allColumnNames = {};
        for i = 1:length(taskFields)
            tableData = S.(taskFields{i}).data;
            if istable(tableData)
                allColumnNames = union(allColumnNames, tableData.Properties.VariableNames, 'stable');
            end
        end

        % Ensure all tables for the current task have the same column names
        for i = 1:length(taskFields)
            tableData = S.(taskFields{i}).data;
            if istable(tableData)
                % Find missing columns for the current table
                missingColumns = setdiff(allColumnNames, tableData.Properties.VariableNames, 'stable');
                
                % Add missing columns and fill them with NaN
                for j = 1:length(missingColumns)
                    tableData.(missingColumns{j}) = NaN(height(tableData), 1);
                end
                
                % Update the table in the new struct
                S_new.(taskFields{i}).data = tableData;
            end
        end
    end
end