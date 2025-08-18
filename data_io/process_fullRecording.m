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
    addParameter(p, 'NASNET_PATH', defaultNET_PATH, @ischar); % only used for plex
    addParameter(p, 'SORTER_PATH', [], @ischar);
    
    % Parse inputs
    parse(p, session_name, varargin{:});
    RAW_PATH    = p.Results.RAW_DATA_PATH;
    OUT_PATH    = p.Results.OUT_DATA_PATH;
    NEV_PATH    = p.Results.NEVUTIL_PATH;
    NET_PATH    = p.Results.NASNET_PATH;
    SORTER_PATH = p.Results.SORTER_PATH;

    addpath(genpath(NEV_PATH));

    % Get the directory of the current script or function
    currentDir = fileparts(mfilename('fullpath'));
    parentDirOneLevelUp = fileparts(currentDir);
    addpath(genpath(parentDirOneLevelUp));
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist(fullfile(OUT_PATH, session_name), 'dir'), mkdir(fullfile(OUT_PATH, session_name)); end

    session_path = fullfile(RAW_PATH, session_name);
    filePattern = fullfile(RAW_PATH, session_name, '*.ns5');

    if isfile(fullfile(session_path,'metadata.json'))
        metadata = loadMetadataJSON(fullfile(session_path,'metadata.json'));
        Tmeta = struct2table(rmfield(metadata, {'sess_name', 'HEeye_VEeye_diode_pupil'}));
    end

    S1 = struct();
    S1.sess_name = metadata.sess_name;

    % Create the search pattern to find files that start with 'filename' and end with '.ns5'
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

    if ismember('neuropixel',metadata.probe_type)
        imec_dirs = dir(fullfile(RAW_PATH, session_name,[session_name, '*_imec*']));
        imec_dirs = arrayfun(@(q) fullfile(q.folder, q.name), imec_dirs, 'uni', 0);
        imec_nums = cellfun(@(q) str2num(q(end)), imec_dirs, 'uni', 0);

        alignCodes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bfv_8_0_9.txt']));
        alignTimes = readmatrix(fullfile(RAW_PATH, session_name, ['catgt_',session_name],[session_name,'_tcat.nidq.bft_8_0_9.txt']));
        alignTimes = alignTimes(alignCodes>0);

        imec_meta = cell(numel(imec_dirs),3);
        for imec = 1:numel(imec_dirs)
            lfp_ap_path = fullfile(imec_dirs{imec}, [session_name, sprintf('_t0.imec%d',imec_nums{imec})]);

            % read in meta data for lfp
            lfp_meta = readMetaFile([lfp_ap_path,'.lf.meta']);
            ap_meta = readMetaFile([lfp_ap_path,'.ap.meta']);

            imec_meta(imec,:) = {(imec), ap_meta, lfp_meta};
        end

        imec_meta = cell2table(imec_meta,'VariableNames', {'probe_index','ap_meta','lfp_meta'});
        S1.metadata = [Tmeta imec_meta];    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Extracting raw data from nev/out datafiles 
    tic

    goodFlag = true;
    for nevnum = 1:length(nevnames) % loop through nev files, in chronological
        nevpath = nevpaths{nevnum};
        this_task = tasks{nevnum};
    
        fprintf('\n---- generating nev_out for %s ----\n', this_task);

        %----- NEUROPIXELS -----%
        if ismember('neuropixel',metadata.probe_type)
            [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', false, 'READ_LFP', false, 'alignPulseEnabled', true);
            startAcquisition = datetime(out_ns5.hdr.timeOrigin, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');

            [dat, ~] = format_datTrials(nev, out_ns5);

            firstSyncPulse = startAcquisition + seconds(dat(1).trialcodes(2,3));
            ripple_pulse_timeStamps = cellfun(@(w) (firstSyncPulse + seconds(w)) - seconds(dat(1).trialcodes(2,3)), cellfun(@(q) q(2,3), {dat.trialcodes}.', 'uni', 0), 'uni', 1);     

            if nevnum == 1
                np_pulse_timeStamps = cellfun(@(w) firstSyncPulse + seconds(w), num2cell(alignTimes - alignTimes(1)), 'uni', 1);
            end
    
            [np_mask, ripple_mask] = match_syncPulses_RipToNP(np_pulse_timeStamps, ripple_pulse_timeStamps);
            fprintf('\n dat has %d rows\n', numel(dat))
            fprintf('np_mask = %d/%d, ripple_mask = %d/%d \n', sum(np_mask), length(np_mask), sum(ripple_mask), length(ripple_mask))  
   
            if isequal(session_name,'kendra_scrappy_0136a_g0') 
                [dat,these_alignTimes,goodFlag] = fix_specificSessions(session_name,np_mask,ripple_mask,alignTimes,dat,goodFlag);
            else
                these_alignTimes = alignTimes(np_mask);
                if sum(np_mask) < length(ripple_mask)
                    dat = dat(ripple_mask);
                    fprintf('\n dat NOW has %d rows', numel(dat))
                end
            end 

            tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);

            % Adding LFP sample to table
            lfp_fs = S1.metadata.lfp_meta(1).imSampRate;
            trial_starts_sec = cellfun(@(q,v) q-(v./1000), num2cell(these_alignTimes), tbl.ALIGN_PULSE(:,1), 'uni', 1);
            trial_ends_sec = trial_starts_sec + (tbl.END_TRIAL./1000);

            trial_starts_lfp_samp = floor(trial_starts_sec * lfp_fs);
            trial_ends_lfp_samp = ceil(trial_ends_sec * lfp_fs);

            tbl.imec_sec = [trial_starts_sec trial_ends_sec];
            tbl.imecLFP_samp = [trial_starts_lfp_samp trial_ends_lfp_samp];


        %----- PLEXON -----%
        elseif ismember('plexon',metadata.probe_type)
            addpath(genpath(NET_PATH));
            
            if exist([nevpath,'.ns2'], 'file') == 2
                [nev, out_ns5, out_ns2] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', true);
                lfp = extract_lfpData(nev,out_ns2,mappings.ripChan_num); 
    
                [dat, ~] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', mappings.ripChan_num);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task, 'LFP', lfp);
            else
                [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'), 'READ_LFP', false);
                [dat, ~] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', mappings.ripChan_num);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
            end

        %----- FHC SINGLE ELECTRODE -----%
        elseif ismember('fhc',metadata.probe_type)
            addpath(genpath(NET_PATH));

            [nev, out_ns5, ~] = extract_nevout(nevpath, 'SPIKE_SORT', true, 'netFolder', fullfile(NET_PATH,'networks'));
            [dat, ~] = format_datTrials(nev, out_ns5, 'NEURAL_CHANNELS', 1);
            tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);

        %----- BEHAVIOR ONLY -----%
        else 
            [nev, out_ns5, ~] = extract_nevout(nevpath);
            if ~isempty(nev)
                [dat, ~] = format_datTrials(nev, out_ns5);
                tbl = convert_smithDat_mayoTbl(dat, 'TASK_NAME', this_task);
            else
                dat = []; tbl = [];
            end
        end
    
        % Convert structures to a cell array of string representations
        if ~isempty(tbl)
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

            S1.(this_task).params = merged_struct;
        end

        S1.(this_task).hdr = out_ns5.hdr;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KILOSORT/NEUROPIXELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if SORTER_PATH
            kilosort_all = []; trlAvg_frs_all = cell(1,numel(imec_dirs)); 
            for imec = 1:numel(imec_dirs)
                %kilosort4_path = fullfile(imec_dirs{imec}, ['kilosort4_', RUN_TYPE]);
                kilosort4_path = fullfile(imec_dirs{imec}, 'sorting', SORTER_PATH);

                if isfolder(kilosort4_path)
                    [spikes_perTrial,kilosort,trlAvg_frs] = parse_KilosortToTbl(tbl,fullfile(kilosort4_path,'sorter_output'),'NP_ALIGN_PULSES',these_alignTimes);
                    kilosort.probe_index = imec;
                    fields = fieldnames(kilosort);
                    fields(strcmp(fields, 'probe_index')) = [];
                    kilosort = orderfields(kilosort, ['probe_index'; fields]);

                    metrics = readtable(fullfile(kilosort4_path,'quality_metrics','cluster_metrics.csv'));
                    kilosort.clusters = [metrics, kilosort.clusters];
                    
                    % kilosort.clusters.sess_name = repmat({metadata.sess_name}, height(kilosort.clusters), 1);
                    kilosort.clusters = movevars(kilosort.clusters,{'sess_name'},'Before','imec');
                    kilosort.clusters.sess_name = categorical(kilosort.clusters.sess_name);

                    trlAvg_frs_all{imec} = trlAvg_frs;
                    kilosort_all = [kilosort_all; kilosort];

                    tbl.(sprintf('spiketimes_%d',imec)) = spikes_perTrial;
                end 
            end

            if nevnum==1
                S1.kilosort = kilosort_all;
            end

            for imec = 1:numel(imec_dirs)
                if ~isempty(trlAvg_frs_all{imec})
                    S1.kilosort(imec).clusters.([this_task, '_Hz']) = trlAvg_frs_all{imec};
                end
            end
 
            %last_alignID = last_alignID + height(tbl);   

            if nevnum==length(nevnames)
                ff = fieldnames(S1);
                S1 = orderfields(S1, ["kilosort"; ff(~strcmp(ff,'kilosort'))]);
            end
        end

        % Remove trials with absolutely no spikes in them
        colnames = {'spiketimes_1', 'spiketimes_2'};
        col_found = colnames(ismember(colnames, tbl.Properties.VariableNames));
        
        if ~isempty(col_found)
            spike_column = tbl.(col_found{1});
            tbl(cellfun(@(q) sum(cellfun(@(w) numel(w), q, 'uni', 1)), tbl.(col_found{1}), 'uni', 1) == 0, :) = [];
        end

        tbl.sess_name = repmat({session_name}, height(tbl), 1);
        tbl = movevars(tbl,{'sess_name'},'Before','trialName');
        tbl.sess_name = categorical(tbl.sess_name);

        S1.(this_task).dat = dat;
        S1.(this_task).tbl = tbl;
    
    end

    % Save the structure S to the specified file
    S = unify_taskTables(S1,taskTypes);

    % if isequal(PROBE_TYPE,'np') && PARSE_KS
    %     S = calculate_metrics_neuropixels(S,trlAvg_frs_all);
    % end

    if ~exist(fullfile(OUT_PATH, session_name, 'tables'), 'dir'), mkdir(fullfile(OUT_PATH, session_name, 'tables')); end 

    if SORTER_PATH
        save(fullfile(OUT_PATH,session_name,'tables',sprintf('%s-%s.mat',session_name,SORTER_PATH)), 'S', '-v7.3');
    else
        save(fullfile(OUT_PATH,session_name,'tables',sprintf('%s.mat',session_name)), 'S', '-v7.3');
    end
    
    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
end
