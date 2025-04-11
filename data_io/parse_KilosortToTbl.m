function [spikes_perTrial,kilosort] = parse_KilosortToTbl(tbl,kilosort4_path,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'kilosort4_path', @ischar);
addParameter(p, 'NP_ALIGN_PULSES', [], @isnumeric);

parse(p, tbl, kilosort4_path, varargin{:});
tbl = p.Results.tbl;
ks_path = p.Results.kilosort4_path;
NP_ALIGN_PULSES = p.Results.NP_ALIGN_PULSES;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in MATLAB table
tic

% Load in kilosort results
kilo_files = dir(fullfile(ks_path, '*.mat'));
for i = 1:length(kilo_files)
    load(fullfile(kilosort4_path,kilo_files(i).name));
end

kilosort.channel_positions = channel_positions;
kilosort.channel_map = channel_map;
kilosort.channel_shanks = channel_shanks;
kilosort.pc_feature_ind = pc_feature_ind;
kilosort.similar_templates = similar_templates;
kilosort.templates = templates;
kilosort.templates_ind = templates_ind;
kilosort.whitening_mat = whitening_mat;
kilosort.whitening_mat_dat = whitening_mat_dat;
kilosort.whitening_mat_inv = whitening_mat_inv;

unique_clusters = double(unique(spike_clusters)+1);

spike_times_sec = double(spike_times)./30000;

spikes_perTrial = cell(height(tbl),1);
for t = 1:height(tbl)
    if mod(t, 100) == 0
        disp(['Trial: ', num2str(t), '/', num2str(height(tbl))]);
    end

    np = NP_ALIGN_PULSES(t);
    rp = tbl.ALIGN_PULSE{t,1};
    et = tbl.END_TRIAL(t);

    spike_units = spike_clusters(spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)));
    spike_times = ((spike_times_sec((spike_times_sec>=(np-(rp./1000)) & spike_times_sec<=(np+((et-rp)./1000)))) - np)*1000) + rp;

    spikes_perTrial{t} = cellfun(@(u) spike_times(spike_units==u), num2cell(unique_clusters), 'uni', 0);
end

end