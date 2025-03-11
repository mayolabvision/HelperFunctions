function process_RippleRecording_KKN(experimenter,monkey,session,varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % experimenter  =  'kendra';
    % monkey        =  'scrappy';
    % session       =  '0097a';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    defaultRAW_PATH  =  '/Volumes/lab_NHPdata';
    defaultOUT_PATH  =  '/Users/kendranoneman/OneDrive/DATA';
    defaultCSV_PATH  =  '/Users/kendranoneman/OneDrive/DATA/RECORDING_INFO.csv';
    defaultNET_PATH  =  '/Users/kendranoneman/Packages/nasnet/networks';
    
    p = inputParser;
    addRequired(p, 'experimenter', @ischar);
    addRequired(p, 'monkey', @ischar);
    addRequired(p, 'session', @ischar);
    addParameter(p, 'SAVE_RAW', True, @islogical); 
    addParameter(p, 'RAW_PATH', defaultRAW_PATH, @ischar); 
    addParameter(p, 'OUT_PATH', defaultOUT_PATH, @ischar); 
    addParameter(p, 'CSV_PATH', defaultCSV_PATH, @ischar); 
    addParameter(p, 'NET_PATH', defaultNET_PATH, @ischar); 
    
    % Parse inputs
    parse(p, experimenter, monkey, session, varargin{:});
    experimenter = p.Results.experimenter;
    monkey = p.Results.monkey;
    session = p.Results.session;
    SAVE_RAW = p.Results.SAVE_RAW;
    RAW_PATH = p.Results.RAW_PATH;
    OUT_PATH = p.Results.OUT_PATH;
    CSV_PATH = p.Results.CSV_PATH;
    NET_PATH   = p.Results.NET_PATH;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    [this_sess,filename]  =  read_recordingNotes(CSV_PATH,experimenter,monkey,session);
    if ~exist(fullfile(OUT_PATH, filename), 'dir'), mkdir(fullfile(OUT_PATH, filename)); end
    
    taskTypes = {'rfmp','purs','mdir','fstm'};
    numTasks = [this_sess.rfmp_num this_sess.purs_num this_sess.mdir_num this_sess.fstm_num];
    
    tasks = [];
    for i=1:length(taskTypes)
        tasks = [tasks, arrayfun(@(x) sprintf('%s%d', taskTypes{i}, x), 1:numTasks(i), 'UniformOutput', false)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Extracting raw data from nev/out datafiles 
    [mappings,probe_specs] = map_channelsNumbersToNames(this_sess.mapFile_name,this_sess.probeID{1},'probeDepths_mm',this_sess.recordDepth_mm{1});
    %mappings.absDepth_mm = mappings.absDepth_mm - (8.3-this_sess.gtHeight_mm{1}{1});
    
    tic
    
    % Make structure to hold all data 
    S1 = struct();
    S1.recording_info = table2struct(this_sess);
    S1.channels = mappings;  
    S1.probe_specs = probe_specs;
    
    start_times = cell(sum(numTasks),1);
    task_num = 1;
    for i = 1:length(taskTypes)
        these_tasks = arrayfun(@(x) sprintf('%s%d', taskTypes{i}, x), 1:numTasks(i), 'UniformOutput', false);
        for f = 1:length(these_tasks)
            this_task = these_tasks{f};
    
            fprintf('\n---- generating nev_out for %s ----\n', this_task);
        
            nevname = sprintf('%s/%s_%s', RAW_PATH, filename, this_task);
    
            if exist(sprintf('%s/%s_%s.ns2', RAW_PATH, filename, this_task), 'file') == 2
                [nev, out_ns5, out_ns2] = extract_nevout(nevname, 'SPIKE_SORT', true, 'netFolder', NET_PATH, 'READ_LFP', true);
                lfp = extract_rawData(nev,out_ns2,mappings.ripChan_num);
    
                tbl = format_dataTable(nev, out_ns5, mappings.ripChan_num, this_task, 'LFP', lfp);
            else
                [nev, out_ns5, ~] = extract_nevout(nevname, 'SPIKE_SORT', true, 'netFolder', NET_PATH, 'READ_LFP', false);
                tbl = format_dataTable(nev, out_ns5, mappings.ripChan_num, this_task);
            end
    
            % Save out_ns5 to raw .bin
            if SAVE_RAW
                if ~contains(this_task,'fstm')
                    bin_path = fullfile(OUT_PATH,filename,[this_task,'.bin']);
                    if exist(bin_path, 'file') == 0
                        fprintf('\n---- writing to bin for %s ----\n', this_task);
                        this_ns5 = out_ns5.data(ismember(out_ns5.hdr.label, string(1:512)),:)';
                        this_ns5 = this_ns5(:,mappings.ripChan_num);
                        
                        fileID = fopen(bin_path, 'wb');
                        fwrite(fileID, this_ns5', 'int16');
                        fclose(fileID);
                    end
                end
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
    
            start_times{task_num} = datetime(out_ns5.hdr.timeOrigin);
            task_num = task_num + 1;
        end
    end
    
    S2 = unify_taskTables(S1,taskTypes);
    
    % Rearrange based on task order
    times = datetime([start_times{:}]);
    [~, idx] = sort(times);
    sorted_tasks = tasks(idx);
    
    S = struct('recording_info', S2.recording_info, 'channels', S2.channels);
    for i = 1:numel(sorted_tasks)
        S.(sorted_tasks{i}) = S2.(sorted_tasks{i});
    end
    
    % Save the structure S to the specified file
    save(fullfile(OUT_PATH,filename,[filename,'.mat']), 'S');
    
    if SAVE_RAW
        % Concatenate .bin files
        full_bin_path = fullfile(OUT_PATH,filename,[filename,'.bin']);
        if exist(full_bin_path, 'file') == 0
            fid_write = fopen(full_bin_path,'w');
            binFiles = dir(fullfile(OUT_PATH, filename, '*.bin'));
            binFiles = binFiles(~strcmp({binFiles.name}, [filename, '.bin']));
        
            for j = 1:length(binFiles)
                fprintf('\n---- concatenating %s to bin ----\n', this_task);
                fid_read = fopen(fullfile(binFiles(j).folder,binFiles(j).name));
                A = fread(fid_read, '*int16');
                fwrite(fid_write, A, 'int16');
                fclose(fid_read);
        
                % Delete individual .bin file after it is written
                delete(fullfile(binFiles(j).folder, binFiles(j).name));
            end
            fclose(fid_write);
        end
    end
    
    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
end