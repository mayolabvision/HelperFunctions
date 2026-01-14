function [mn,err,yu,yl] = get_errorBars_overTime(X, varargin)
% ERRORBAR_BOUNDS Mean and error bounds (SEM or STD)
%
% [mn,err,yu,yl] = errorbar_bounds(X)
% [mn,err,yu,yl] = errorbar_bounds(X,'ErrorType','std','OmitNaN',true)
%
% Inputs:
%   X          - NxM array (observations x time/conditions)
%
% Name-Value:
%   'ErrorType' - 'sem' (default) or 'std'
%   'OmitNaN'   - true or false (default = false)
%
% Outputs:
%   mn  - mean across rows
%   err - SEM or STD
%   yu  - upper bound (mn + err)
%   yl  - lower bound (mn - err)

p = inputParser;
addParameter(p,'ErrorType','sem',@(s) any(strcmpi(s,{'sem','std'})));
addParameter(p,'OmitNaN',false,@islogical);
parse(p,varargin{:});

errType = lower(p.Results.ErrorType);
omitnan = p.Results.OmitNaN;

% --- mean ---
if omitnan
    mn = mean(X,1,'omitnan');
else
    mn = mean(X,1);
end

% --- error ---
switch errType
    case 'sem'
        if omitnan
            n  = sum(~isnan(X),1);
            err = std(X,'omitnan') ./ sqrt(n);
        else
            err = std(X) ./ sqrt(size(X,1));
        end
    case 'std'
        if omitnan
            err = std(X,'omitnan');
        else
            err = std(X);
        end
end

yu = mn + err;
yl = mn - err;

end