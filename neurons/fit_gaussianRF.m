function [params, gof] = fit_gaussianRF(xvals, yvals, meanFR)
% fit_gaussianRF — Fit a 2D Gaussian to a receptive field firing-rate map.
%
% INPUTS:
%   xvals   : 1xNx vector of x (horizontal) grid positions
%   yvals   : 1xNy vector of y (vertical) grid positions
%   meanFR  : Ny x Nx matrix of mean firing rate (Hz) at each grid position,
%             rows = y, cols = x (matches format_tableToRFMap's convention)
%
% OUTPUT:
%   params  : [amp, x0, y0, sigX, sigY, offset] of the fitted Gaussian
%               FR(x,y) = amp * exp(-((x-x0)^2/(2*sigX^2) + (y-y0)^2/(2*sigY^2))) + offset
%   gof     : struct with .rsquare, .sse, .ssTot, and .hitBound. hitBound is
%             true if the optimizer parked at a constraint instead of
%             converging to an interior solution -- amp ~ 0 (no real
%             response), or sigX/sigY pinned to their (generous) upper cap
%             (not localized). x0/y0 and sigma's LOWER bound are excluded:
%             a real RF's center can legitimately sit just outside the
%             tested stimulus range (an eccentric RF the display couldn't
%             fully capture), and a real RF can be narrower than one grid
%             cell (the grid just can't resolve anything finer) -- neither
%             means the fit failed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X, Y] = meshgrid(xvals, yvals); % Ny x Nx, matches meanFR orientation
XY = cat(3, X, Y);

offset0 = min(meanFR(:));
amp0 = max(meanFR(:)) - offset0;
[~, peakIdx] = max(meanFR(:));
[peakRow, peakCol] = ind2sub(size(meanFR), peakIdx);

dx = abs(mean(diff(xvals))); dy = abs(mean(diff(yvals)));
params0 = [amp0, xvals(peakCol), yvals(peakRow), range(xvals)/4, range(yvals)/4, offset0];
lb = [0, min(xvals)-range(xvals), min(yvals)-range(yvals), dx/2, dy/2, -Inf];
ub = [Inf, max(xvals)+range(xvals), max(yvals)+range(yvals), range(xvals)*2, range(yvals)*2, Inf];

gaussFun = @(p, XY) p(1) * exp(-(((XY(:,:,1)-p(2)).^2)/(2*p(4)^2) + ((XY(:,:,2)-p(3)).^2)/(2*p(5)^2))) + p(6);

opts = optimoptions('lsqcurvefit', 'Display', 'off');
params = lsqcurvefit(gaussFun, params0, XY, meanFR, lb, ub, opts);

fitted = gaussFun(params, XY);
ssRes = sum((meanFR(:)-fitted(:)).^2);
ssTot = sum((meanFR(:)-mean(meanFR(:))).^2);
gof.sse = ssRes;
gof.ssTot = ssTot;
gof.rsquare = 1 - ssRes/ssTot;

tol = 1e-3;
nearBound = @(v, b) abs(v-b) < tol*max(abs(b), 1);
gof.hitBound = nearBound(params(1), lb(1)) ...      % amp ~ 0, no real response
    || nearBound(params(4), ub(4)) ...              % sigX pinned to the upper cap -- not localized
    || nearBound(params(5), ub(5));                 % sigY pinned to the upper cap -- not localized

end
