function convert_ns5ToBin(session_name, varargin)
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

    %defaultRAW_PATH  =  '/Volumes/lab_NHPdata';
    defaultRAW_PATH  =  '/Volumes/home/DATA';

    p = inputParser;
    addRequired(p, 'session_name', @ischar);
    addParameter(p, 'RAW_DATA_PATH', defaultRAW_PATH, @ischar); 
    addParameter(p, 'CHANNELS', [], @isnumeric); 

    parse(p, session_name, varargin{:});
    RAW_PATH = p.Results.RAW_DATA_PATH;
    CHANNELS = p.Results.CHANNELS; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filePattern = fullfile(RAW_PATH, session_name, '*.ns5');
    raw_files = dir(filePattern);
    raw_filenames = {raw_files.name}.';
    nevnames = cellfun(@(q) q(1:end-4), raw_filenames, 'uni', 0);
    raw_filepaths = arrayfun(@(x) fullfile(x.folder, x.name), raw_files, 'UniformOutput', false);

    recording_times = cellfun(@(l) l.hdr.timeOrigin, cellfun(@(q) read_nsx(q,'readdata',false), raw_filepaths, 'uni', 0), 'uni', 0);
    [~,idx] = sort(recording_times);
    nevnames = nevnames(idx);
    nevpaths = raw_filepaths(idx);

    % Define possible task keywords
    task_keywords = {'rfmp', 'rfMapping', 'purs', 'pursuit', 'mdir', 'dirmem', 'fstm'};
    
    % Initialize cell array for tasks
    tasks = cell(size(nevnames));
    
    % Loop through each nevname
    for i = 1:numel(nevnames)
        name = nevnames{i};
        found = false;
        for j = 1:numel(task_keywords)
            pattern = [task_keywords{j}, '\w*'];  % keyword followed by letters/numbers
            match = regexp(name, pattern, 'match', 'once');
            if ~isempty(match)
                tasks{i} = match;
                found = true;
                break;
            end
        end
    end

    tic

    full_bin_path = fullfile(RAW_PATH, session_name, [session_name, '.bin']);
    if exist(full_bin_path, 'file') == 2
        delete(full_bin_path);
    end

    for nevnum = 1:length(nevnames)
        nevpath = nevpaths{nevnum};
        this_task = tasks{nevnum};

        fprintf('\n---- generating nev_out for %s ----\n', this_task);
        
        [~, out_ns5, ~] = extract_nevout(nevpath);

        if ~contains(this_task,'fstm')
            this_ns5 = out_ns5.data(ismember(out_ns5.hdr.label, string(1:512)),:);
            if ~isempty(this_ns5)
                fprintf('\n---- writing to bin for %s ----\n', this_task);
    
                if ~isempty(CHANNELS)
                    this_ns5 = this_ns5(CHANNELS,:);
                end

                fid_write = fopen(full_bin_path, 'a'); % Open file in append mode ('a')
                fwrite(fid_write, this_ns5, 'int16');
    
                fclose(fid_write);  
            else
                fprintf('\n---- no raw signal for %s ----\n', this_task);
            end
        end
    end
    

end