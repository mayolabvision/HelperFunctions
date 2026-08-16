function [unitnames,snrs,exp_clean] = findConsistentUnits_fromStruct(exp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    exp_clean = exp;

    if isfield(exp_clean.info, 'channels')
        channels       =  exp_clean.info.channels; % names of all channels
    else
        units  =  cellfun(@(x) fieldnames(x), {exp_clean.dataMaestroPlx.units}.', 'uni', 0);
        channels   =  sort(unique(vertcat(units{:})));
    end
    channels = cellfun(@(y) strrep(y,'unit',''), channels, 'uni', 0);
    if isfield(exp_clean.info, 'SNRs')
        snrs  =  exp_clean.info.SNRs; % SNR for each channel
    else
        snrs  =  nan(1,length(channels));
    end
     
    % throw out trials that are missing too many units
    %all_units  =  cellfun(@(x) fieldnames(x), {exp_clean.dataMaestroPlx.units}.', 'uni', 0);
    %exp_clean.dataMaestroPlx(cellfun(@(q) size(q,1), all_units) <= min(cellfun(@(q) size(q,1), all_units))) = [];

    all_units  =  cellfun(@(x) fieldnames(x), {exp_clean.dataMaestroPlx.units}.', 'uni', 0);
    [B,BG]         =  groupcounts(vertcat(all_units{:}));
    [~,ia]         =  setdiff(channels,cellfun(@(y) strrep(y,'unit',''), BG(B==max(B)), 'uni', 0));
    channels(ia)   =  []; snrs(ia) = [];
    
    [unitnames,I]  =  sort(channels); snrs = snrs(I);
    unitnames      =  cellfun(@(z) strcat('unit',z), unitnames, 'uni', 0)';
end