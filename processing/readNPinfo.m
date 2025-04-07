function [trialinfo, spikesbyclust] = readNPinfo(session)
%% NP event and spike pre-processing
% read in and merge alignment pulses from NP files, spike times/clusters,
% and nev event codes

% all filename pulls assume we are in lab or neuropixels directory

% relevant files:
% NP alignment info (.txt), output from catGT
% NP spikes: spike_times.npy,
% spikes_clusters.npy (if no manual sorting done), or
% cluster_group.tsv (if manual sorting was done)
% .nev file for trial code info

% output: spikesbytrial: spike times per trial, organized by cluster,
% and event codes per trial
% trialinfo:
% col1: trial dig codes aligned to pulse
% col2: trial direction
% spikes by cluster: spike times aligned to pulse and organized by cluster
% and trial

% session = 's244';
% block = '0003';

%% get NP alignment info
% extract info from txt files: NP align times and values (0 or 1, time in sec)

% works whether NP alignment pulses are in separate files or one

files = dir(fullfile(strcat('neuropixels/walter_data/*', session, '*/*bfv', '*.txt')));
fileIDs = {files.name};
alignNP_vals = [];
for  i = 1:length(fileIDs)
    name = fileIDs{i};
    alignNP_vals = vertcat(alignNP_vals, fscanf(fopen(name), '%f'));
end

files = dir(fullfile(strcat('neuropixels/walter_data/*', session, '*/*bft', '*.txt')));
fileIDs = {files.name};
alignNP_times = [];
for i = 1:length(fileIDs)
    name = fileIDs{i};
    alignNP_times = vertcat(alignNP_times, fscanf(fopen(name), '%f'));
end

alignNP = [alignNP_vals, alignNP_times];

%% find breaks in blocks - ie where one block ended and another began
% find difference between consecutive elements
diffs = abs(diff((alignNP(:, 2))));

% Find indices where the difference is greater than 10s (this is where
% block changes occur)
idx = find(diffs > 10);

idx = [idx, idx+1]; %col 1 is end of block, col2 is beginning of next block
blockstarts = alignNP(idx(:,2),2); %times a new behav block begins

%% process all nevs for alignment and events
files = dir(fullfile(strcat('neuropixels/walter_data/*', session, '*/*.nev')));
filenames = {files.name};
alldat = {};

for i = 1:length(filenames)
    alldat{i} = nev2dat(char(filenames(i)));
end
% correct trials only
% datcorrect = dat([dat.result] == 150);

% align first trial of recording: align to code 1 - make this time 0 of
% both nev events and NP file...

%% get spikes
% get clusters/times from KS output
% create spikes array:
% col1: cluster assignment; col2: times (in s)

% convert spike_times and spike_clusters from npy to mat
files = dir(fullfile(strcat('neuropixels/walter_data/*', session, '*/kilosort4/spike_clusters.npy')));
fileIDs = files.name;
spikes = readNPY(fileIDs);

files = dir(fullfile(strcat('neuropixels/walter_data/*', session, '*/kilosort4/spike_times.npy')));
fileIDs = files.name;
spikes(:,2) = readNPY(fileIDs);

% time in samples, 30kHz samp rate
% convert to double and get time in s
spikes = double(spikes);
spikes(:,2) = spikes(:,2)/30000;



%% option1: align times of nev codes and NP alignment times, transfer those times to spike times
%
% codes = {dat.trialcodes};
% codes = vertcat(codes{:});

%% alt: divide spikes into trials and blocks based on align times

%{
%spikesbytrial(numtrials).clusters = [];
spikesbytrial(numtrials).spikes = [];
%spikesbytrial(numtrials).trialcodes = [];


% TODO: for now, not dealing with last behav. trial; there is no "end" code
% in NP aligns
for i = 1:numtrials-1
    idx = spikes(:,2)>trialstarts(i, 2) & spikes(:,2)<trialstarts(i+1, 2);
    % spikesbytrial(i).clusters = spikes(idx, 1);
    % spikesbytrial(i).times = spikes(idx, 2);
    spikesbytrial(i).trialcodes = dat(i).trialcodes;
    %spikesbytrial(i).spikes = spikes(idx, :);
    
    trialspks = spikes(idx, :);
    trialclusts = unique(trialspks(:,1));
    spikesbytrial(i).spikes = cell(length(trialclusts), 2);

    % try to organize spikes field by cluster...
    for j = 1:length(trialclusts)
        clust = trialclusts(j);
        clusttimes = trialspks(trialspks(:,1) == clust, 2);
        spikesbytrial(i).spikes{j,1} = clust;
        spikesbytrial(i).spikes{j, 2} = clusttimes;
    end

end
%}

%% trialcodes aligned to pulse,
% realign trialcode and spike times per trial where pulse-on = time 0
% each 1 in alignNP is start of a trial

numblocks = length(alldat);
trialstarts = alignNP(alignNP(:,1) == 1, :);

tcodesaligned = cell(1, numblocks);

for iblock = 1:numblocks
    dat = alldat{iblock};
    numtrials = length(dat);
    t = cell(numtrials,1);
    for itrial = 1:numtrials
        % make time 0 = pulse-on time (will be same for spike times)
        tcodes = dat(itrial).trialcodes;
        talign = tcodes(tcodes(:,2)==0, 3);
        tcodes(:,3)=tcodes(:,3)-talign(1);
        t{itrial} = tcodes;
    end
    tcodesaligned{1, iblock} = t;
end


%% spikesbycluster
allclusters = unique(spikes(:,1));
spikesbyclust(length(allclusters)).cluster = [];

for iclust = 1:length(allclusters)
    clust = allclusters(iclust);
    spikesbyclust(iclust).cluster = clust;
    clustspks = spikes(spikes(:,1)==clust, 2); %spike times for cluster
    % clust spike times by trial
    for iblock = 1:numblocks
        dat = alldat{iblock};
        name = strcat('spikes_block', num2str(iblock));
        spikesbyclust(iclust).(name) = [];
        if iblock == 1
            numtrials = length(dat);
            trial1 = 1;
        else
            trial1 = numtrials + 1;
            numtrials = numtrials + length(dat);
        end

        for itrial = trial1:numtrials
            blocktrial = itrial+1-trial1;
            tcodes = tcodesaligned{1,iblock}{blocktrial};

            trialend = tcodes(tcodes(:,2)==255, 3);
            idx = clustspks>trialstarts(itrial, 2) & clustspks<trialstarts(itrial, 2)+trialend;
            trialspks = clustspks(idx);
            talignNP = trialstarts(itrial,2);
            spikesbyclust(iclust).(name){blocktrial, 1} = trialspks-talignNP;

            % spike times, make 0 = pulse-on time (same for dig codes)

        end
    end
end

%% extract direction info (conditions) for each trial

trialinfo = cell(1,numblocks);

for iblock = 1:numblocks
    dat = alldat{iblock};
    directions = ones(length(dat),1);
    delays = ones(length(dat),1);
    trialinfo{iblock}.dirs = [];
    trialinfo{iblock}.delays = [];
    trialinfo{iblock}.tcodes = tcodesaligned{iblock};
    for itrial = 1:length(dat)
        conds = dat(itrial).text;
        angle = strfind(conds, "angle");
        fix = strfind(conds, "fix");
        trialdir = str2double(conds(angle+6:fix-2));
        trialdelay = str2num(conds(fix+12:end));
        directions(itrial) = trialdir;
        delays(itrial) = trialdelay;
    end
    trialinfo{iblock}.dirs = directions;
    trialinfo{iblock}.delays = delays;
end


%% to extract spike times for all trials for one cluster of interest:

% spikesbyclust([spikesbyclust.cluster]==297);


end

