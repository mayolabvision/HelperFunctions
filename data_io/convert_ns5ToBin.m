function convert_ns5ToBin(ns5, save_path, varargin)
% convert_ns5ToBin: Extracts raw neural signals from NS5 files or matrices and writes to a .bin file.
%
% Supports input as:
%   - A single NS5 file path (char)
%   - A struct containing NS5 data
%   - A cell array of NS5 file paths (concatenates in order)
%
% Inputs:
%   - ns5: (char, struct, or cell array of char) NS5 data input.
%   - save_path: (char) Path to save the binary file.
%
% Optional Parameters:
%   - 'CHANNELS': (numeric, default = all channels)
%   - 'CHUNK_SIZE': (numeric, default = 1e6) Rows processed at a time.
%
% Example:
%   convert_ns5ToBin({'file1.ns5', 'file2.ns5'}, 'output.bin');

p = inputParser;
addRequired(p, 'ns5', @(x) (ischar(x) || isstruct(x) || iscell(x)));
addRequired(p, 'save_path', @(x) ischar(x) || isstring(x));
addParameter(p, 'CHANNELS', [], @isnumeric); 

% Parse inputs
parse(p, ns5, save_path, varargin{:});
CHANNELS = p.Results.CHANNELS;

% Open file for writing
fid_write = fopen(save_path, 'w');

for j = 1:length(ns5)
    fid_read = fopen
end


    for i = 1:numel(ns5)
        fprintf('\n---- binning raw data for task %d/%d ----\n', i, numel(ns5));
        write_ns5_to_bin(ns5{i}, fileID, CHANNELS, CHUNK_SIZE);
    end
else
    write_ns5_to_bin(ns5, fileID, CHANNELS, CHUNK_SIZE);
end

% Close file
fclose(fileID);

end

function write_ns5_to_bin(ns5, fileID, CHANNELS, CHUNK_SIZE)
    % Reads NS5 file and writes it to an open binary file.
    
    if isstruct(ns5)
        this_ns5 = ns5.data;
    elseif ischar(ns5)
        [~, out_ns5, ~] = extract_nevout(ns5, 'SPIKE_SORT', false, 'READ_LFP', false);
        this_ns5 = out_ns5.data;
    else
        error('Unsupported data type for ns5 input.');
    end
    
    % Ensure data is in correct orientation
    if size(this_ns5,1) < size(this_ns5,2)
        this_ns5 = this_ns5';
    end
    
    this_ns5 = this_ns5(:, ismember(out_ns5.hdr.label, string(1:512)));
    if ~isempty(CHANNELS)
        this_ns5 = this_ns5(:, CHANNELS);
    end
    
    

    % Convert data to int16
    this_ns5 = int16(this_ns5);
    
    % Write data in chunks with a progress bar
    numRows = size(this_ns5, 1);
    numChunks = ceil(numRows / CHUNK_SIZE);
    % h = waitbar(0, 'Writing data...');
    
    for chunkIdx = 1:numChunks
        startIdx = (chunkIdx - 1) * CHUNK_SIZE + 1;
        endIdx = min(chunkIdx * CHUNK_SIZE, numRows);
        
        % Write chunk to file
        fwrite(fileID, this_ns5(startIdx:endIdx, :), 'int16');
        
    end
    

end