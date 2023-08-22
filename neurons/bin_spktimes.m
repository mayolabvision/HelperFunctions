function spks_binned = bin_spktimes(spktimes,startTime,endTime,binSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

edges = startTime : binSize : endTime;

allCounts = zeros(1, length(edges) - 1);
[counts, values] = histcounts(spktimes, edges);
spks_binned = allCounts + counts;

end