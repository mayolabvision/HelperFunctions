function d = circDiffDeg(a, b, signed)
% circDiffDeg  Circular difference between angles in degrees
%
% d = circDiffDeg(a, b)
% d = circDiffDeg(a, b, signed)
%
% Inputs:
%   a, b   : scalars or arrays of angles (degrees), same size
%   signed : (optional) true → signed diff in [-180, 180]
%            false (default) → unsigned diff in [0, 180]
%
% Output:
%   d      : circular difference (same size as a and b)

if nargin < 3
    signed = false;
end

d = mod(a - b + 180, 360) - 180;

if ~signed
    d = abs(d);
end
end