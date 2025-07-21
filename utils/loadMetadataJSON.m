function metadata = loadMetadataJSON(filename)
    % Load a JSON file into a struct if it exists
    if isfile(filename)
        fid = fopen(filename, 'r');
        raw = fread(fid, inf, 'char=>char')';
        fclose(fid);
        metadata = jsondecode(raw);
    else
        warning('File %s not found.', filename);
        metadata = struct();  % return empty struct if file doesn't exist
    end
end