function [hasRF, params, rsq, pValue] = analyzeRFDirection(meanFRs, xvals, yvals, binIdx, direction, ALPHA, MIN_RSQ, MAX_AMP)
% analyzeRFDirection — Test one direction ('exc' or 'inh') of a receptive
% field at one time bin: fits a 2D Gaussian (to the map itself for 'exc', to
% the inverted map for 'inh'), then tests significance via a nested-model
% F-test comparing that fit (6 params) against a flat/no-RF null (1 param:
% the grand mean), over the grid's cell means. This has far more power to
% detect real, smooth, modest-firing-rate RFs than a per-position omnibus
% ANOVA over raw per-flash spike counts, which dilutes the signal across
% many largely-noisy single-flash comparisons.
%
% INPUTS:
%   meanFRs   : Ny x Nx x N_BINS array of mean firing rate (Hz)
%   xvals     : 1xNx vector of x (horizontal) grid positions
%   yvals     : 1xNy vector of y (vertical) grid positions
%   binIdx    : which of the N_BINS time bins to test
%   direction : 'exc' or 'inh'
%   ALPHA     : significance threshold for the F-test
%   MIN_RSQ   : minimum Gaussian fit R^2 for a candidate RF to count
%   MAX_AMP   : max plausible |amplitude| (Hz) for a candidate RF to count
%               (an implausibly large amplitude signals an outlier-driven
%               spurious fit rather than real sustained firing)
%
% OUTPUTS:
%   hasRF   : true if the F-test is significant, R^2 > MIN_RSQ, |amplitude|
%             <= MAX_AMP, and the fit didn't park at an optimizer bound
%             instead of converging (fit_gaussianRF.m's gof.hitBound)
%   params  : [amp, x0, y0, sigX, sigY, offset] of the fitted Gaussian
%   rsq     : the fit's R^2
%   pValue  : the F-test p-value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peakSlice = meanFRs(:,:,binIdx);

if strcmp(direction, 'exc')
    [params, gof] = fit_gaussianRF(xvals, yvals, peakSlice);
else
    [paramsRaw, gof] = fit_gaussianRF(xvals, yvals, -peakSlice); % invert so the dip becomes a peak to fit
    params = paramsRaw;
    params(1) = -paramsRaw(1); % report as negative amplitude = depth of suppression
    params(6) = -paramsRaw(6);
end
rsq = gof.rsquare;

n = numel(peakSlice);
pFit = 6; pNull = 1;
F = ((gof.ssTot - gof.sse) / (pFit - pNull)) / (gof.sse / (n - pFit));
pValue = 1 - fcdf(F, pFit - pNull, n - pFit);

hasRF = pValue < ALPHA && gof.rsquare > MIN_RSQ && abs(params(1)) <= MAX_AMP && ~gof.hitBound;

end
