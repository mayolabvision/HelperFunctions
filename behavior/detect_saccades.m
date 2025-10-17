function saccades = detect_saccades(eyeVel,varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    p = inputParser;
    addRequired(p, 'eyeVel', @isnumeric);
    addParameter(p, 'VEL_THRESH', 30, @isnumeric)
    addParameter(p, 'ACC_THRESH', 1000, @isnumeric)

    parse(p, eyeVel, varargin{:});

    eyeVel = p.Results.eyeVel;
    VEL_THRESH = p.Results.VEL_THRESH;
    ACC_THRESH = p.Results.ACC_THRESH;
    
    if size(eyeVel,2) > size(eyeVel,1)
        vel = eyeVel;
    else
        vel = eyeVel';
    end

    if size(vel,1) == 1
        rVel = vel;
    else
        [~,rVel] = cart2pol(vel(1,:),vel(2,:));
    end

    % calculate acceleration 
    dt = 1 / 1000;                     % time per sample (s)
    offset = round(0.010 * 1000);      % number of samples in 10 ms
    
    % Use central difference: (v[t+10ms] - v[t-10ms]) / 20ms
    rAcc = (circshift(rVel(1, :), -offset) - circshift(rVel(1, :), offset)) / (2 * offset * dt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1. Identify candidate saccade regions
    saccMask = (abs(rVel) > VEL_THRESH) | (abs(rAcc) > ACC_THRESH);
    
    % Find starts and ends of 1-blocks
    d = diff([0, saccMask, 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    
    runLengths = ends - starts + 1;
    validBlocks = (runLengths >= 6) & starts > 10 & starts < (length(rVel)-20);
    
    starts = starts(validBlocks);
    saccOnsets = nan(size(starts));
    for i = 1:numel(starts)
        saccOnsets(i) = find(rAcc(1:starts(i)) < ACC_THRESH, 1, 'last') + 1;
    end
    saccOnsets = unique(saccOnsets);

    saccOffsets = nan(size(saccOnsets));
    for i = 1:numel(saccOnsets)
        if i ~= numel(saccOnsets)
            [~,m] = max(rVel(saccOnsets(i):saccOnsets(i+1)));
        else
            [~,m] = max(rVel(saccOnsets(i):end));
        end

        if ~isempty(find(diff(sign(rAcc(saccOnsets(i)+m:end))) ~= 0, 1, 'first'))
            saccOffsets(i) = saccOnsets(i) + find(diff(sign(rAcc(saccOnsets(i)+m:end))) ~= 0, 1, 'first');
        end
    end

    if isnan(saccOffsets(i))
        saccOnsets = saccOnsets(1:end-1);
        saccOffsets = saccOffsets(1:end-1);
    end

    blocks = [saccOnsets; saccOffsets]';

    if ~isempty(blocks)
        mergeGap = 50;  % max gap to consider merging
    
        merged = [];
        currentStart = blocks(1,1);
        currentEnd   = blocks(1,2);
        
        for i = 2:size(blocks,1)
            nextStart = blocks(i,1);
            nextEnd   = blocks(i,2);
            
            if nextStart - currentEnd <= mergeGap
                currentEnd = max(currentEnd, nextEnd);
            else
                merged = [merged; currentStart, currentEnd]; 
                currentStart = nextStart;
                currentEnd   = nextEnd;
            end
        end
        
        % Add the last block
        saccades = [merged; currentStart, currentEnd];
    else
        saccades = blocks;
    end

end