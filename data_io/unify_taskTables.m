function S_new = unify_taskTables(S, tasknames)
    % This function unifies task-specific tables within a structure by ensuring 
    % that all tables related to the same task have the same set of columns. 
    % It takes a structure S containing multiple fields, each holding task data 
    % in the form of tables, and a cell array of tasknames. For each task, it 
    % finds all tables associated with that task, identifies the union of all 
    % column names across those tables, and ensures that each table contains 
    % all the columns, filling missing ones with NaN values. The function returns 
    % a new structure, S_new, with updated tables where all task-related tables 
    % have consistent columns.

    % INPUTS:
    % S - A structure containing fields, each corresponding to a task with 
    %     table data.
    % tasknames - A cell array of task names, where each task name corresponds 
    %             to a set of related fields in the structure S.

    % OUTPUTS:
    % S_new - A new structure with updated tables where each task's related 
    %         tables have the same columns.

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
            tableData = S.(taskFields{i}).tbl;
            if istable(tableData)
                allColumnNames = union(allColumnNames, tableData.Properties.VariableNames, 'stable');
            end
        end

        % Ensure all tables for the current task have the same column names
        for i = 1:length(taskFields)
            tableData = S.(taskFields{i}).tbl;
            if istable(tableData)
                % Find missing columns for the current table
                missingColumns = setdiff(allColumnNames, tableData.Properties.VariableNames, 'stable');
                
                % Add missing columns and fill them with NaN
                for j = 1:length(missingColumns)
                    tableData.(missingColumns{j}) = NaN(height(tableData), 1);
                end
                
                % Update the table in the new struct
                S_new.(taskFields{i}).tbl = tableData;
            end
        end
    end
end
