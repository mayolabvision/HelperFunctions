function [mn,err,yu,yl] = sem_errorbar(X, varargin)

% default
errType = 'sem';

if ~isempty(varargin)
    errType = lower(varargin{1});
end

% mean
mn = mean(X,1);

% sample size
n = size(X,1);

% std
sd = std(X,0,1);

switch errType
    case 'sem'
        err = sd ./ sqrt(n);

    case 'std'
        err = sd;

    case {'ci','ci95','95ci'}
        err = 1.96 * (sd ./ sqrt(n));

    otherwise
        error('errType must be ''sem'', ''std'', or ''ci95''')
end

yu = mn + err;
yl = mn - err;

end