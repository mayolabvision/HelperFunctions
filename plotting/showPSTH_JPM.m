function [psth, newRange, sampling] = showPSTH_JPM(spikeTimes,bounds,sigma, col)
%function psth = showPSTH(spikeTimes,bounds,sigma)
%
% Code stolen from the Smith Lab utilities by JPM on July 7 2014

% spikeTimes is a list of spikeTimes (in s)
%  if spikeTimes is a cell array, all values inside will be pooled
%
%sigma, if included, is the standard deviation of the gaussian smoothing window
%
% col = color, default is black

if nargin < 4
    col = 'k'; %color
end

  sampling = .0001;
  
if iscell(spikeTimes)
    spikeTimes = reshape(spikeTimes, numel(spikeTimes), 1); % dump all spikes times into a single vertical vector
    trials = length(spikeTimes);
    
    for i = 1:length(spikeTimes)
        spikeTimes{i} = spikeTimes{i}(:);
    end
    
    newRange = [bounds(1)-5*sigma*sampling bounds(2)+5*sigma*sampling];
    
    spikeTimes = cell2mat(spikeTimes);
    spikeTimes = spikeTimes(spikeTimes >= bounds(1)-5*sigma*sampling & spikeTimes <= bounds(2)+5*sigma*sampling);
    
    psth = hist(spikeTimes,newRange(1):sampling:newRange(2));
else
  trials = size(spikeTimes,1);  
  
  psth = sum(spikeTimes,1);
  newRange = sampling*[1 length(psth)];
  bounds = [newRange(1) newRange(2)];
end

psth = gaussSmooth(psth,sigma);
psth = psth/sampling/trials;

samps = length(bounds(1)+sampling:sampling:bounds(2));

if nargout == 0
    plot(newRange(1):sampling:newRange(2),psth,col)
    xlim([bounds(1) bounds(2)])
end

toCut = (length(psth)-samps+1)/2;
psth = psth(toCut:end-toCut);