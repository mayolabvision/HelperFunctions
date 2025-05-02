function process_fullRecording(session_name,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % session_name = 'kendra_scrappy_0136a_g0' 
    % session_name is the name of a datafolder, which contains ripple data
    % (.ns5, .nev, etc...) and folders of spikeglx imec probe data, catgt align pulses 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default paths to add to the MATLAB path    
    defaultRAW_PATH  =  '/Volumes/lab_NHPdata';
    defaultOUT_PATH  =  '/Volumes/home/DATA';
    defaultNEV_PATH  =  '/Users/kendranoneman/Packages/nevutils';
    defaultNET_PATH  =  '/Users/kendranoneman/Packages/nasnet';
    
    p = inputParser;
    addRequired(p, 'session_name', @ischar);
    addParameter(p, 'RAW_DATA_PATH', defaultRAW_PATH, @ischar); 
    addParameter(p, 'OUT_DATA_PATH', defaultOUT_PATH, @ischar); 
    addParameter(p, 'NEVUTIL_PATH', defaultNEV_PATH, @ischar);
    addParameter(p, 'PROBE_TYPE', [], @ischar); % np, plex
    addParameter(p, 'NASNET_PATH', defaultNET_PATH, @ischar); % only used for plex
    addParameter(p, 'PARSE_KILOSORT', false, @islogical);
    addParameter(p, 'RUN_TYPE', 'unleashed', @ischar);
    
    % Parse inputs
    parse(p, session_name, varargin{:});
    RAW_PATH = p.Results.RAW_DATA_PATH;
    OUT_PATH = p.Results.OUT_DATA_PATH;
    NEV_PATH = p.Results.NEVUTIL_PATH;
    PROBE_TYPE = p.Results.PROBE_TYPE;
    NET_PATH = p.Results.NASNET_PATH;
    PARSE_KS = p.Results.PARSE_KILOSORT;
    RUN_TYPE = p.Results.RUN_TYPE;

    addpath(genpath(NEV_PATH));

    % Get the directory of the current script or function
    currentDir = fileparts(mfilename('fullpath'));
    parentDirOneLevelUp = fileparts(currentDir);
    addpath(genpath(parentDirOneLevelUp));
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist(fullfile(OUT_PATH, session_name), 'dir'), mkdir(fullfile(OUT_PATH, session_name)); end

    % Create the search pattern to find files that start with 'filename' and end with '.ns5'
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
        if ~found
            tasks{i} = 'unknown';
        end
    end

    taskTypes = unique(cellfun(@(q) regexp(q, '[a-zA-Z]+', 'match', 'once'), tasks, 'uni', 0));

    disp(tasks)
    if isequal(PROBE_TYPE,'np')
        imec_dirs = dir(fullfile(OUT_PATH, session_name,[session_name, '*_imec*']));
        imec_dirs = arrayfun(@(q) fullfile(q.folder, q.name), imec_dirs, 'uni', 0);
        imec_nums = cellfun(@(q) str2num(q(end)), imec_dirs, 'uni', 0);

        alignCodes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bfv_8_0_9.txt']));
        alignTimes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bft_8_0_9.txt']));
        alignTimes = alignTimes(alignCodes>0);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Extracting raw data from nev/out datafiles 
    tic
    
    S1 = struct();

    lastFileEnd = 0; last_alignID = 0; goodFlag = true;
    for nevnum = 1:length(nevnames)
        nevpath = nevpaths{nevnum};
        this_task = tasks{nevnum};
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIPPLE BEHAVIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n---- generating nev_out for %s ----\n', this_task);

        if isequal(PROBE_TYPE, 'np')
            [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false, 'alignPulseEnabled', true);
            startAcquisition = datetime(out_ns5.hdr.timeOrigin, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');

            [dat, nsEnd] = format_datTrials(nev, out_ns5);

            firstSyncPulse = startAcquisition + seconds(dat(1).trialcodes(2,3));
            ripple_pulse_timeStamps = cellfun(@(w) (firstSyncPulse + seconds(w)) - seconds(dat(1).trialcodes(2,3)), cellfun(@(q) q(2,3), {dat.trialcodes}.', 'uni', 0), 'uni', 1);     

            if nevnum == 1
                np_pulse_timeStamps = cellfun(@(w) firstSyncPulse + seconds(w), num2cell(alignTimes - alignTimes(1)), 'uni', 1);
            end
    
            [np_mask, ripple_mask] = match_syncPulses_RipToNP(np_pulse_timeStamps, ripple_pulse_timeStamps);
    
    
            if goodFlag % only used for kendra_scrappy_0136a_g0
                if sum(np_mask) >= numel(ripple_mask)
                    these_alignTimes = alignTimes(np_mask);
                    goodFlag = true;
        
                elseif sum(np_mask) < numel(ripple_mask) % only used for kendra_scrappy_0136a_g0
                    first_block_start = find(np_mask, 1, 'first');
                    first_block_end = first_block_start + find(~np_mask(first_block_start:end), 1, 'first') - 2;
                    first_zero_index = first_block_end + 1;
     
                    good_alignTimes1 = alignTimes(first_block_start:first_block_end);
                    remaining_alignTimes = alignTimes(first_zero_index:end);
                    good_alignTimes2 = remaining_alignTimes(1:696);
        
                    these_alignTimes = [good_alignTimes1; good_alignTimes2];
                    dat(772:869) = [];
        
                    goodFlag = false;
                end
            else
                these_alignTimes = alignTimes(end-313:end);
            end

            tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);

        elseif isequal(PROBE_TYPE, 'plex')
            if exist([nevpath,'.ns2'], 'file') == 2
                [nev, out_ns5, out_ns2] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', true);
                lfp = extract_lfpData(nev,out_ns2,mappings.ripChan_num); 
    
                [dat, nsEnd] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', mappings.ripChan_num);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task, 'LFP', lfp);
            else
                [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', false);
                [dat, nsEnd] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', mappings.ripChan_num);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
            end

        end

        % add to ns5_samps
        % if ~contains(this_task,'fstm')
        %     tbl.ns5_samps = cellfun(@(q) q+lastFileEnd, num2cell(tbl.ns5_samps,2), 'uni', 0);
        %     lastFileEnd = lastFileEnd + nsEnd;
        % else
        %     tbl.ns5_samps = [];
        % end
    
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

            kilosort_all = []; trlAvg_frs_all = cell(1,numel(imec_dirs));
            for imec = 1:numel(imec_dirs)
                kilosort4_path = fullfile(imec_dirs{imec}, ['kilosort4_', RUN_TYPE]);

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

            for imec = 1:numel(imec_dirs)
                S1.kilosort(imec).clusters.([this_task, '_Hz']) = trlAvg_frs_all{imec};
            end
            
            last_alignID = last_alignID + height(tbl);   

            if nevnum==length(nevnames)
                ff = fieldnames(S1);
                S1 = orderfields(S1, ["kilosort"; ff(~strcmp(ff,'kilosort'))]);
            end
        end

        S1.(this_task).data = tbl;
    
    end

    % Save the structure S to the specified file
    S = unify_taskTables(S1,taskTypes);
    
    save(fullfile(OUT_PATH,session_name,sprintf('%s_%s.mat',session_name,RUN_TYPE)), 'S', '-v7.3');
    
    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
end
