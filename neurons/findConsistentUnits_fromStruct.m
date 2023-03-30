function [unitnames,snrs] = findConsistentUnits_fromStruct(exp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    channels       =  exp.info.channels; % names of all channels
    snrs           =  exp.info.SNRs; % SNR for each channel 
    
    all_units      =  cellfun(@(x) fieldnames(x), {exp.dataMaestroPlx.units}.', 'uni', 0);
    [B,BG]         =  groupcounts(vertcat(all_units{:}));
    [~,ia]         =  setdiff(channels,cellfun(@(y) y(end-3:end), BG(B==max(B)), 'uni', 0));
    channels(ia)   =  []; snrs(ia) = [];
    
    [unitnames,I]  =  sort(channels); snrs = snrs(I);
    unitnames      =  cellfun(@(z) strcat('unit',z), unitnames, 'uni', 0)';
end