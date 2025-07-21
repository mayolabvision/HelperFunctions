function T = metadataStructToTable(metadata)
    % Converts a struct with scalar/mixed-length fields to a table
    % Repeats singleton fields to match the length of the longest field

    fields = fieldnames(metadata);
    maxLength = 1;

    % First determine the number of rows (based on longest field)
    for i = 1:numel(fields)
        val = metadata.(fields{i});
        if iscell(val) || isnumeric(val)
            len = size(val, 1);  % assumes column vectors or cell arrays
        else
            len = 1;
        end
        maxLength = max(maxLength, len);
    end

    % Now build a cell array row-by-row
    cellOut = cell(maxLength, numel(fields));
    for i = 1:numel(fields)
        val = metadata.(fields{i});

        if iscell(val)
            if isscalar(val)
                val = repmat(val, maxLength, 1);
            elseif size(val,1) ~= maxLength
                error('Field %s has incompatible size.', fields{i});
            end
            cellOut(:,i) = val;

        elseif isnumeric(val)
            if isscalar(val)
                cellOut(:,i) = repmat({val}, maxLength, 1);
            elseif isvector(val) && size(val,1) == maxLength
                cellOut(:,i) = num2cell(val);
            else
                error('Field %s has incompatible size.', fields{i});
            end

        elseif ischar(val) || isstring(val)
            cellOut(:,i) = repmat({val}, maxLength, 1);

        else
            cellOut(:,i) = repmat({val}, maxLength, 1);
        end
    end

    % Convert to table
    T = cell2table(cellOut, 'VariableNames', fields);
end