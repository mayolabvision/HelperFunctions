function [spikes_perTrial,kilosort,trlAvg_frs] = parse_KilosortToTbl(tbl,kilosort4_path,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'kilosort4_path', @ischar);
addParameter(p, 'NP_ALIGN_PULSES', [], @isnumeric);

parse(p, tbl, kilosort4_path, varargin{:});
tbl = p.Results.tbl;
ks_path = p.Results.kilosort4_path;
NP_ALIGN_PULSES = p.Results.NP_ALIGN_PULSES;

tokens = regexp(ks_path, 'imec(\d+)', 'tokens');
if isempty(tokens)
    error('Could not find "imec" followed by a number in ks_path');
end
imec_num = str2double(tokens{1}{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in MATLAB table
tic
kilo_files1 = dir(fullfile(ks_path, '*.mat'));
for i = 1:length(kilo_files1)
    kpath = fullfile(kilosort4_path, kilo_files1(i).name);
    [~, kname, ~] = fileparts(kpath);
    
    loaded = load(kpath);  % Load into struct
    varNames = fieldnames(loaded);  % Get the name(s) of loaded variables

    loaded_var = loaded.(varNames{1});
    if isa(loaded_var, 'single') || isa(loaded_var, 'int32') || isa(loaded_var, 'int64')
        loaded_var = double(loaded_var);
    end
    
    % Assign the first variable in the .mat file to the struct
    kilosort.(kname) = loaded_var;
end

% .tsv files
kilo_files2 = dir(fullfile(ks_path, '*.tsv'));
[~, sort_idx] = sort({kilo_files2.name}.');
kilo_files2 = kilo_files2(sort_idx);

for ii = 1:length(kilo_files2)
    cfpath = fullfile(kilo_files2(ii).folder, kilo_files2(ii).name);
    [~, cfname, ~] = fileparts(kilo_files2(ii).name);
    if ii==1
        clusters = readtable(cfpath, 'FileType', 'text', 'Delimiter', '\t');
    else
        cc = readtable(cfpath, 'FileType', 'text', 'Delimiter', '\t');
        cc2 = join(clusters, cc, 'Keys', 'cluster_id');
        clusters = cc2;
    end
end
clusters.KSLabel_clusters = categorical(clusters.KSLabel_clusters);

clusters.probe_index = repmat(imec_num+1, height(clusters), 1);
kilosort.clusters = movevars(clusters, 'probe_index', 'Before', 1);

unique_clusters = double(unique(kilosort.spike_clusters)+1);
spike_times_sec = double(kilosort.spike_times)./30000;

spikes_perTrial = cell(height(tbl),1);
for t = 1:height(tbl)
    if mod(t, 100) == 0
        disp(['Trial: ', num2str(t), '/', num2str(height(tbl))]);
    end

    np = NP_ALIGN_PULSES(t);
    rp = tbl.ALIGN_PULSE{t,1};
    et =  tbl.END_TRIAL(t);

    spike_units = kilosort.spike_clusters(spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)));
    spike_times = ((spike_times_sec((spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)))) - np)*1000) + rp;

    spikes_perTrial{t} = cellfun(@(u) spike_times(spike_units==u), num2cell(unique_clusters), 'uni', 0);
end

spike_counts = cellfun(@(w) cellfun(@(q) numel(q), w, 'uni', 1), spikes_perTrial, 'uni', 0);
unit_frs = cellfun(@(w,v) w./v, spike_counts, num2cell((tbl.END_TRIAL - tbl.START_TRIAL)./1000), 'uni', 0);
unit_frs = unit_frs(cellfun(@sum, unit_frs, 'uni', 1)>0);
trlAvg_frs = mean(vertcat(unit_frs{:}),1,'omitnan')';

end
