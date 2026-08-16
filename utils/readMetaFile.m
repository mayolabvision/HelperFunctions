function meta = readMetaFile(filename)
    % Reads a SpikeGLX .meta file and returns a struct of key-value pairs
    meta = struct();

    fid = fopen(filename, 'rt');
    if fid == -1
        error('Failed to open meta file: %s', filename);
    end

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || startsWith(line, '#')
            continue;  % Skip empty or comment lines
        end

        tokens = split(line, '=');
        if numel(tokens) == 2
            key = strtrim(tokens{1});
            valueStr = strtrim(tokens{2});

            % Attempt numeric conversion
            valueNum = str2double(valueStr);
            if ~isnan(valueNum)
                value = valueNum;
            else
                % Try parsing comma-separated list of numbers
                if contains(valueStr, ',')
                    items = split(valueStr, ',');
                    nums = str2double(items);
                    if all(~isnan(nums))
                        value = nums;
                    else
                        value = valueStr;  % Keep as string
                    end
                else
                    value = valueStr;  % Keep as string
                end
            end

            meta.(matlab.lang.makeValidName(key)) = value;
        end
    end

    fclose(fid);
end