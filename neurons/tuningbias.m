function [dirsel, dirpref, orisel, oripref, quadsel, quadpref] = tuningbias(theta, rho, bl)
%function [dirsel, dirpref, orisel, oripref, quadsel, quadpref] = tuningbias(theta, rho, bl)
%
% uses vector averaging to calculate tuning stats
% theta is angle in degrees, rho is response magnitude, bl is
% spont activity baseline (optional - 0 if not passed)
%
% Calculates a tuning bias vector (Leventhal et al., 1995; O’Keefe and
% Movshon, 1998), similar to the vector strength calculation introduced
% by Levick and Thibos (1982).

% angle conversion functions
%r2d = inline('((x./pi).*180)');  % radians to degrees
%d2r = inline('(pi.*(x./180))');  % degrees to radians

if nargin < 3
        bl = 0;
end;

% make into column vectors
rho = rho(:);
theta = theta(:);

% subtract baseline from responses
rho = rho - bl;
% sum response magnitudes
srho = sum(abs(rho));

%%%%%%%%%%%%%%%%%%%%%%%
% quad tuning
%%%%%%%%%%%%%%%%%%%%%%%

% vector sum using complex exponentials
% theta multiplied by 4 to "stretch" around circle 4x to
% get quad selectivity (irrespective of orientation)
sel = sum(rho.*exp(i.*4.*deg2rad(theta)));

% selectivity = length of vector sum / sum of vector lengths
quadsel = abs(sel)./rho;

% preferred ori is direction of vector sum (/2 since we stretched it)
quadpref = mod(rad2deg(angle(sel))/2,90);

%%%%%%%%%%%%%%%%%%%%%%%
% orientation tuning
%%%%%%%%%%%%%%%%%%%%%%%

% vector sum using complex exponentials
% theta multiplied by 2 to "stretch" around circle twice to
% get orientation selectivity (irrespective of direction)
sel = sum(rho.*exp(i.*2.*deg2rad(theta)));

% selectivity = length of vector sum / sum of vector lengths
orisel = abs(sel)./srho;

% preferred ori is direction of vector sum (/2 since we stretched it)
oripref = mod(rad2deg(angle(sel))/2,180);

%%%%%%%%%%%%%%%%%%%%%%%
% direction tuning
%%%%%%%%%%%%%%%%%%%%%%%

% vector sum using complex exponentials
sel = sum(rho.*exp(i.*deg2rad(theta)));

% selectivity = length of vector sum / sum of vector lengths
dirsel = abs(sel)./srho;

% preferred ori is direction of vector sum
dirpref = mod(rad2deg(angle(sel)),360);

