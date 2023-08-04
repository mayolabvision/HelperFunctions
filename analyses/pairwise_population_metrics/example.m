clear; clc; close all;
addpath(genpath('/Users/kendranoneman/Projects/mayo/pairwise_population_metrics'))

%% load data

% the two spike count matrices (counts.attend_in, counts.attend_out) are of
% the shape n_neurons, n_trials
% load('example.mat');
counts = struct; 
counts.saccadeTask = table2array(stbl.saSpks_zscore{stbl.Session=='06A' & stbl.AlignmentType=='stim' & stbl.TimeWindowStart==350});
counts.pursuitTask = table2array(stbl.puSpks_zscore{stbl.Session=='06A' & stbl.AlignmentType=='stim' & stbl.TimeWindowStart==350});

% because pairwise and population metrics can depend on n_neurons and n_trials,
% we recommend that, if possible, one use the same neurons and equalize the 
% number of trials between conditions one is comparing
[n_neurons,n_trials] = size(counts.saccadeTask);
fprintf('number of neurons: %d\n',n_neurons)
fprintf('number of trials: %d\n\n',n_trials)

%% compute pairwise metrics

% compute the mean and s.d. of the rsc distribution
[rsc_mean_sa,rsc_sd_sa] = compute_pairwise_metrics(counts.saccadeTask);
[rsc_mean_pu,rsc_sd_pu] = compute_pairwise_metrics(counts.pursuitTask);

% compare the two conditions
fprintf('Rsc mean:\n')
fprintf('   saccade: %.3f, pursuit: %.3f\n\n',rsc_mean_sa,rsc_mean_pu)

fprintf('Rsc standard deviation:\n')
fprintf('   saccade: %.3f, pursuit: %.3f\n\n',rsc_sd_sa,rsc_sd_pu)

%% compute population metrics

% fit factor analysis & compute population metrics:
%   1) perform cross-validation to choose dimensionalities among those 
%      specified by 'zDimList'
%   2) compute population metrics & shared eigenspectrum
% the function 'compute_pop_metrics' performs both of these steps

zDimList = 5; %0:15;

[ls_sa,psv_sa,d_sa,espec_sa] = ...
    compute_population_metrics(counts.saccadeTask,zDimList);
[ls_pu,psv_pu,d_pu,espec_pu] = ...
    compute_population_metrics(counts.pursuitTask,zDimList);

% compare the two conditions
fprintf('Loading similarity of the strongest dimension:\n')
fprintf('   saccade: %.3f, pursuit: %.3f\n\n',ls_sa(1),ls_pu(1))

fprintf('%% shared variance:\n')
fprintf('   saccade: %.3f, pursuit: %.3f\n\n',psv_sa,psv_pu)

fprintf('Dimensionality (d_shared):\n')
fprintf('   saccade: %d, pursuit: %d\n\n',d_sa,d_pu)
