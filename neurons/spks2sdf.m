function sdf = spks2sdf(spks, time, sigma)
% spks  : vector of spike times
% time  : time axis
% sigma : Gaussian width

sdf = arrayfun(@(mu) ...
    exp(-0.5 * ((time - mu) / sigma).^2) ./ (sigma * sqrt(2*pi)), ...
    spks, 'uni', 0);

end