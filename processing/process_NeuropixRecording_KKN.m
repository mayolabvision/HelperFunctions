function process_NeuropixRecording_KKN(experimenter,monkey,session,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % experimenter  =  'kendra';
    % monkey        =  'scrappy';
    % session       =  '0097a';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default paths to add to the MATLAB path    
    defaultRAW_PATH  =  '/Volumes/lab_NHPdata';
    defaultOUT_PATH  =  '/Volumes/home/DATA';
    defaultCSV_PATH  =  '/Volumes/home/RECORDING_INFO.csv';
    defaultNPY_PATH  =  '/Users/kendranoneman/Packages/npy-matlab';
    defaultNEV_PATH  =  '/Users/kendranoneman/Packages/nevutils';
    
    p = inputParser;
    addRequired(p, 'experimenter', @ischar);
    addRequired(p, 'monkey', @ischar);
    addRequired(p, 'session', @ischar);
    % addParameter(p, 'SAVE_RAW', true, @islogical); 
    addParameter(p, 'RAW_DATA_PATH', defaultRAW_PATH, @ischar); 
    addParameter(p, 'OUT_DATA_PATH', defaultOUT_PATH, @ischar); 
    addParameter(p, 'RECD_CSV_PATH', defaultCSV_PATH, @ischar); 
    addParameter(p, 'NPY_MAT_PATH', defaultNPY_PATH, @ischar);
    addParameter(p, 'NEVUTIL_PATH', defaultNEV_PATH, @ischar); 
    
    % Parse inputs
    parse(p, experimenter, monkey, session, varargin{:});
    experimenter = p.Results.experimenter;
    monkey = p.Results.monkey;
    session = p.Results.session;
    % SAVE_RAW = p.Results.SAVE_RAW;
    RAW_PATH = p.Results.RAW_DATA_PATH;
    OUT_PATH = p.Results.OUT_DATA_PATH;
    CSV_PATH = p.Results.RECD_CSV_PATH;
    NPY_PATH = p.Results.NPY_MAT_PATH;
    NEV_PATH = p.Results.NEVUTIL_PATH;

    addpath(genpath(NPY_PATH));
    addpath(genpath(NEV_PATH));

    % Get the directory of the current script or function
    currentDir = fileparts(mfilename('fullpath'));
    parentDirOneLevelUp = fileparts(currentDir);
    addpath(genpath(parentDirOneLevelUp));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [this_sess,filename]  =  read_recordingNotes(CSV_PATH,experimenter,monkey,session);
    if ~exist(fullfile(OUT_PATH, filename), 'dir'), mkdir(fullfile(OUT_PATH, filename)); end

    % Create the search pattern to find files that start with 'filename' and end with '.ns5'
    filePattern = fullfile(RAW_PATH, [filename,'_g0'],[filename, '*.ns5']);
    raw_files = dir(filePattern);
    raw_filenames = {raw_files.name}.';
    nevnames = cellfun(@(q) q(1:end-4), raw_filenames, 'uni', 0);
    raw_filepaths = arrayfun(@(x) fullfile(x.folder, x.name), raw_files, 'UniformOutput', false);

    recording_times = cellfun(@(l) l.hdr.timeOrigin, cellfun(@(q) read_nsx(q,'readdata',false), raw_filepaths, 'uni', 0), 'uni', 0);
    [~,idx] = sort(recording_times);
    nevnames = nevnames(idx);
    nevpaths = raw_filepaths(idx);

    tasks = cellfun(@(q) q{4}, cellfun(@(x) split(x, '_'), nevnames, 'uni', 0), 'uni', 0);
    taskTypes = unique(cellfun(@(q) regexp(q, '[a-zA-Z]+', 'match', 'once'), tasks, 'uni', 0));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Extracting raw data from nev/out datafiles 
    % [mappings,probe_specs] = map_channelsNumbersToNames(this_sess.mapFile_name,this_sess.probeID{1},'probeDepths_mm',this_sess.recordDepth_mm{1});
    % mappings.absDepth_mm = mappings.absDepth_mm - (8.3-this_sess.gtHeight_mm{1}{1});
    
    tic
    
    % Make structure to hold all data 
    S1 = struct();
    S1.recording_info = table2struct(this_sess);
    % S1.channels = mappings;  
    % S1.probe_specs = probe_specs;

    % % Check if the file already exists, if so delete
    % if SAVE_RAW
    %     full_bin_path = fullfile(OUT_PATH, filename, [filename, '.bin']);
    %     if exist(full_bin_path, 'file') == 2
    %         delete(full_bin_path);
    %     end
    % end

    lastFileEnd = 0;
    for nevnum = 1:length(nevnames)
        nevpath = nevpaths{nevnum};
        this_task = tasks{nevnum};
    
        fprintf('\n---- generating nev_out for %s ----\n', this_task);
    
        [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false, 'alignPulseEnabled', true);
        [tbl,nsEnd] = format_dataTable(nev, out_ns5, 'CONVERT_TO_TABLE', true, 'TASK_NAME', this_task);

        % add to ns5_samps
        if ~contains(this_task,'fstm')
            tbl.ns5_samps = cellfun(@(q) q+lastFileEnd, num2cell(tbl.ns5_samps,2), 'uni', 0);
            lastFileEnd = lastFileEnd + nsEnd;
        else
            tbl.ns5_samps = [];
        end
    
        % Convert structures to a cell array of string representations
        all_params = {tbl.params.block}.';
        structStrings = cellfun(@(x) jsonencode(x), all_params, 'UniformOutput', false);
        [~, uniqueIdx] = unique(structStrings, 'stable');
        unique_structs = all_params(uniqueIdx);
        merged_struct = struct();
        fieldNames = fieldnames(unique_structs{1});
        for ii = 1:numel(fieldNames)
            field = fieldNames{ii};
            merged_struct.(field) = cellfun(@(s) s.(field), unique_structs, 'UniformOutput', false);
            if all(cellfun(@isnumeric, merged_struct.(field)))
                merged_struct.(field) = cell2mat(merged_struct.(field));
            end
        end

        S1.(this_task).hdr = out_ns5.hdr;
        S1.(this_task).params = merged_struct;
        S1.(this_task).data = tbl;
    
    end

    
    % Save the structure S to the specified file
    S = unify_taskTables(S1,taskTypes);
    save(fullfile(OUT_PATH,filename,[filename,'.mat']), 'S');
    
    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
end
