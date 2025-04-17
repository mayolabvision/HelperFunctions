function process_NeuropixRecording_KKN(session_name,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % experimenter  =  'kendra';
    % monkey        =  'scrappy';
    % session       =  '0097a';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default paths to add to the MATLAB path    
    defaultRAW_PATH  =  '/Volumes/lab_NHPdata';
    defaultOUT_PATH  =  '/Volumes/home/DATA';
    defaultCSV_PATH  =  '/Volumes/home/DATA/RECORDING_INFO.csv';
    defaultNEV_PATH  =  '/Users/kendranoneman/Packages/nevutils';
    
    p = inputParser;
    addRequired(p, 'session_name', @ischar);
    addParameter(p, 'RAW_DATA_PATH', defaultRAW_PATH, @ischar); 
    addParameter(p, 'OUT_DATA_PATH', defaultOUT_PATH, @ischar); 
    addParameter(p, 'RECD_CSV_PATH', defaultCSV_PATH, @ischar); 
    addParameter(p, 'NEVUTIL_PATH', defaultNEV_PATH, @ischar); 
    addParameter(p, 'PARSE_KILOSORT', true, @islogical)
    addParameter(p, 'DRIFT_CORRECTED', false, @islogical)
    
    % Parse inputs
    parse(p, session_name, varargin{:});
    RAW_PATH = p.Results.RAW_DATA_PATH;
    OUT_PATH = p.Results.OUT_DATA_PATH;
    CSV_PATH = p.Results.RECD_CSV_PATH;
    NEV_PATH = p.Results.NEVUTIL_PATH;
    PARSE_KS = p.Results.PARSE_KILOSORT;
    DRIFT = p.Results.DRIFT_CORRECTED;

    addpath(genpath(NEV_PATH));

    % Get the directory of the current script or function
    currentDir = fileparts(mfilename('fullpath'));
    parentDirOneLevelUp = fileparts(currentDir);
    addpath(genpath(parentDirOneLevelUp));
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fparts = split(session_name, '_');
    experimenter = fparts{1}; monkey = fparts{2}; session = fparts{3};

    [this_sess,~]  =  read_recordingNotes(CSV_PATH,experimenter,monkey,session);
    if ~exist(fullfile(OUT_PATH, session_name), 'dir'), mkdir(fullfile(OUT_PATH, session_name)); end

    % Create the search pattern to find files that start with 'filename' and end with '.ns5'
    filePattern = fullfile(RAW_PATH, session_name, [session_name(1:end-3), '*.ns5']);
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
    
    if PARSE_KS
        imec_dirs = dir(fullfile(OUT_PATH, session_name,[session_name, '*_imec*']));
        imec_dirs = arrayfun(@(q) fullfile(q.folder, q.name), imec_dirs, 'uni', 0);

        imec_nums = cellfun(@(q) str2num(q(end)), imec_dirs, 'uni', 0);
    end

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

    alignCodes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bfv_8_0_9.txt']));
    alignTimes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bft_8_0_9.txt']));
    alignTimes = alignTimes(alignCodes==2);

    lastFileEnd = 0; last_alignID = 0;
    for nevnum = 1:length(nevnames)
        nevpath = nevpaths{nevnum};
        this_task = tasks{nevnum};
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIPPLE BEHAVIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n---- generating nev_out for %s ----\n', this_task);

        [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false, 'alignPulseEnabled', true);
        [dat, nsEnd] = format_datTrials(nev, out_ns5);
        tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEUROPIXELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if PARSE_KS
            % SYNC PULSE
            these_alignTimes = alignTimes(last_alignID+1:last_alignID+height(tbl));

            kilosort_all = []; trlAvg_frs_all = cell(1,numel(imec_dirs));
            for imec = 1:numel(imec_dirs)
                if DRIFT
                    kilosort4_path = fullfile(imec_dirs{imec},'corrected','kilosort4');
                else
                    kilosort4_path = fullfile(imec_dirs{imec},'kilosort4');
                end
                if isfolder(kilosort4_path)
                    [spikes_perTrial,kilosort,trlAvg_frs] = parse_KilosortToTbl(tbl,kilosort4_path,'NP_ALIGN_PULSES',these_alignTimes);
                    trlAvg_frs_all{imec} = trlAvg_frs;
                    kilosort_all = [kilosort_all; kilosort];

                    tbl.(sprintf('spiketimes_imec%d',imec_nums{imec})) = spikes_perTrial;
                end
            end

            if nevnum==1
                S1.kilosort = kilosort_all;
            end

            %for imec = 1:numel(imec_dirs)
            %    S1.kilosort(imec).clusters.([this_task, '_Hz']) = trlAvg_frs_all{imec};
            %end
            
            last_alignID = last_alignID + height(tbl);   
        end

        S1.(this_task).data = tbl;
    
    end

    %ff = fieldnames(S1);
    %S1 = orderfields(S1, [ff(1:find(strcmp(ff,'recording_info'))); "kilosort"; ff(~strcmp(ff,'kilosort') & ~strcmp(ff,'recording_info'))]);

    % Save the structure S to the specified file
    S = unify_taskTables(S1,taskTypes);
    S.rfmp1.data = S.rfmp1.data(~cellfun(@(q) any(isnan(q)), S.rfmp1.data.STIM_OFF, 'uni', 1),:);

    if DRIFT
        save(fullfile(OUT_PATH,session_name,[session_name,'_corrected.mat']), 'S', '-v7.3');
    else
        save(fullfile(OUT_PATH,session_name,[session_name,'.mat']), 'S', '-v7.3');
    end 
    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
end
