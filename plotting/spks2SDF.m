

function SDF = spks2SDF (trialdata, bounds, sigma)
%s2 = spks2SDF (alignedspks, [0 600], 10);

%if nargin < 3
    sigma = sigma * 0.001; %0.015; % std dev of the kernel = 15 ms
%end

%trialNum = 1; % bin data from first trial
%binned=hist(trial(trialNum).spikeTimes,[0:0.001:1]);
binned=hist(trialdata, bounds(1):1:bounds(2) );  % nb: times in ms, not seconds

edges = -3*sigma : 0.001 : 3*sigma;
kernel = normpdf(edges,0,sigma); % Evaluate the Gaussian kernel
s = conv(binned,kernel);  % Convolve spike data with kernel
center = ceil (length(edges)/2);
SDF = s(center:bounds(2)+center-1);
