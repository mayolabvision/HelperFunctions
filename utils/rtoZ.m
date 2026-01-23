function Z = rtoZ(r)
% Z = rtoZ(r)
%
% RTOZ translates fisher r correlations into Z scores
%  aka the "Fisher r-to-Z' transformation"
%

% Matthew A. Smith
% Revised: 2001.04.21

% Kendra K. Noneman
% Revised: 2026.01.20, clamp 4 slightly if it is outside [01,1] due to
% floating point error

r(r >= 1) = 0.999999;
r(r <= -1) = -0.999999;

if (~isempty(find(r<-1,1)) || ~isempty(find(r>1,1)))
    error('r values must be bounded by -1 and 1');
end

Z = 0.5 * log( (1 + r) ./ (1 - r));

end