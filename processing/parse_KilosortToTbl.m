function tbl = parse_KilosortToTbl(experimenter,monkey,session,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultOUT_PATH  =  '/Users/kendranoneman/OneDrive/DATA';

p = inputParser;
addRequired(p, 'experimenter', @ischar);
addRequired(p, 'monkey', @ischar);
addRequired(p, 'session', @ischar);
addParameter(p, 'OUT_DATA_PATH', defaultOUT_PATH, @ischar);

parse(p, experimenter, monkey, session, varargin{:});
experimenter = p.Results.experimenter;
monkey = p.Results.monkey;
session = p.Results.session;
OUT_PATH = p.Results.OUT_DATA_PATH;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in MATLAB table
tic
data_path = fullfile(OUT_PATH,sprintf('%s_%s_%s',experimenter,monkey,session));
load(fullfile(data_path,sprintf('%s_%s_%s.mat',experimenter,monkey,session)));

% Load in kilosort results
kilo_files = dir(fullfile(fullfile(data_path,'kilosort4'), '*.mat'));
for i = 1:length(kilo_files)
    load(fullfile(data_path,'kilosort4',kilo_files(i).name));
end

S.channels.channel_positions = channel_positions;
S.kilosort.pc_feature_ind = pc_feature_ind;
S.kilosort.similar_templates = similar_templates;
S.kilosort.templates = templates;
S.kilosort.templates_ind = templates_ind;
S.kilosort.whitening_mat = whitening_mat;
S.kilosort.whitening_mat_dat = whitening_mat_dat;
S.kilosort.whitening_mat_inv = whitening_mat_inv;

unique_clusters = double(unique(spike_clusters)+1);

S_fnames = fieldnames(S);
for j = 1:length(S_fnames)
    if isfield(S.(S_fnames{j}),'data') && ~contains(S_fnames{j},'fstm')
        fprintf('\n---- parsing kilosort for %s ----\n', S_fnames{j});
        this_task = S.(S_fnames{j}).data;
        this_task.spiketimes = [];
        this_task.diode = [];
        this_task.net_labels = [];
        this_task.eyePos_raw = [];

        %[spiketimes_ms,spike_amps,pc_feats,spike_detect_temps,spike_pos,spike_temps] = deal(cell(height(this_task),length(unique_clusters)));
        spiketimes_ms= cell(height(this_task),length(unique_clusters));
        for c = 1:length(unique_clusters)
            spike_inds =  cellfun(@(q) spike_clusters==unique_clusters(c) & (spike_times>=q(1) & spike_times<=q(2)), this_task.ns5_samps, 'uni', 0);
            spiketimes_ms(:,c) = cellfun(@(q,w) (double(spike_times(q)-w(1))./30000)*1000, spike_inds, this_task.ns5_samps, 'uni', 0);
            % spike_amps(:,c) =  cellfun(@(q) double(amplitudes(q)), spike_inds, 'uni', 0);
            % pc_feats(:,c) =  cellfun(@(q) double(pc_features(q,:,:)), spike_inds, 'uni', 0);
            % spike_detect_temps(:,c) =  cellfun(@(q) double(spike_detection_templates(q)), spike_inds, 'uni', 0);
            % spike_pos(:,c) =  cellfun(@(q) double(spike_positions(q,:)), spike_inds, 'uni', 0);
            % spike_temps(:,c) =  cellfun(@(q) double(spike_templates(q)), spike_inds, 'uni', 0);
        end

        this_task.ks_spiketimes = spiketimes_ms;
        % this_task.ks_spikeamps = spike_amps;
        % this_task.ks_pcFeatures = pc_feats;
        % this_task.ks_spikeDetectTemplates = spike_detect_temps;
        % this_task.ks_spikePositions = spike_pos;
        % this_task.ks_spikeTemplates = spike_temps; 

        S.(S_fnames{j}).data = this_task;      
    end
end

save(fullfile(data_path,sprintf('%s_%s_%s_KS.mat',experimenter,monkey,session)), 'S');

tc = toc;
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

end