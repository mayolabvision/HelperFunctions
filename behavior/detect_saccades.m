function saccades = detect_saccades(eyeVel, varargin)
    % DETECT_SACCADES  Detect saccade onset and offset times from eye velocity.
    %
    %   saccades = detect_saccades(eyeVel)
    %   saccades = detect_saccades(eyeVel, 'Fs', 1000, 'VEL_THRESH', 30, 'ACC_THRESH', 1000)
    %
    %   Inputs:
    %     eyeVel      - Eye velocity trace. Either a 1xN vector (speed) or a
    %                   2xN / Nx2 matrix of [Vx; Vy] components, in deg/s.
    %
    %   Optional name-value parameters:
    %     Fs          - Sampling rate in Hz (default: 1000)
    %     VEL_THRESH  - Velocity threshold in deg/s (default: 30)
    %     ACC_THRESH  - Acceleration threshold in deg/s^2 (default: 1000)
    %
    %   Output:
    %     saccades    - Mx2 matrix of [onset, offset] sample indices, one row
    %                   per detected saccade.

    p = inputParser;
    addRequired(p, 'eyeVel', @isnumeric);
    addParameter(p, 'VEL_THRESH', 30, @isnumeric)
    addParameter(p, 'ACC_THRESH', 1000, @isnumeric)
    addParameter(p, 'Fs', 1000, @isnumeric)

    parse(p, eyeVel, varargin{:});

    eyeVel     = p.Results.eyeVel;
    VEL_THRESH = p.Results.VEL_THRESH;
    ACC_THRESH = p.Results.ACC_THRESH;
    Fs         = p.Results.Fs;

    % --- Ensure velocity is stored row-major (components x samples) ----------
    if size(eyeVel, 2) > size(eyeVel, 1)
        vel = eyeVel;
    else
        vel = eyeVel';
    end

    % Collapse 2D (Vx, Vy) to radial speed; pass through if already 1D
    if size(vel, 1) == 1
        rVel = vel;
    else
        [~, rVel] = cart2pol(vel(1,:), vel(2,:));
    end

    % --- Acceleration via central difference (v[t+10ms] - v[t-10ms]) / 20ms -
    dt     = 1 / Fs;                    % time per sample (s)
    offset = round(0.010 * Fs);         % number of samples in 10 ms

    rAcc = (circshift(rVel(1,:), -offset) - circshift(rVel(1,:), offset)) / (2 * offset * dt);

    % --- 1. Build a binary mask of candidate saccade samples -----------------
    % A sample is flagged if it exceeds either threshold.
    saccMask = (abs(rVel) > VEL_THRESH) | (abs(rAcc) > ACC_THRESH);

    % Find contiguous runs of flagged samples
    d      = diff([0, saccMask, 0]);
    starts = find(d ==  1);
    ends   = find(d == -1) - 1;

    % Discard runs that are too short or too close to the signal boundaries
    runLengths  = ends - starts + 1;
    validBlocks = (runLengths >= 6) & starts > 10 & starts < (length(rVel) - 20);
    starts = starts(validBlocks);

    % --- 2. Refine onsets: walk back to where acceleration first exceeded threshold
    saccOnsets = nan(size(starts));
    for i = 1:numel(starts)
        saccOnsets(i) = find(rAcc(1:starts(i)) < ACC_THRESH, 1, 'last') + 1;
    end
    saccOnsets = unique(saccOnsets);

    % --- 3. Find offsets: locate the zero-crossing of acceleration after peak velocity
    saccOffsets = nan(size(saccOnsets));
    for i = 1:numel(saccOnsets)
        % Search for peak velocity between this onset and the next
        if i ~= numel(saccOnsets)
            [~, m] = max(rVel(saccOnsets(i):saccOnsets(i+1)));
        else
            [~, m] = max(rVel(saccOnsets(i):end));
        end

        % Offset is the first sign change of acceleration after the velocity peak
        zeroCross = find(diff(sign(rAcc(saccOnsets(i)+m:end))) ~= 0, 1, 'first');
        if ~isempty(zeroCross)
            saccOffsets(i) = saccOnsets(i) + zeroCross;
        end
    end

    % Drop the last saccade if its offset could not be determined
    if isnan(saccOffsets(i))
        saccOnsets  = saccOnsets(1:end-1);
        saccOffsets = saccOffsets(1:end-1);
    end

    blocks = [saccOnsets; saccOffsets]';

    % --- 4. Merge saccades separated by a short gap --------------------------
    if ~isempty(blocks)
        mergeGap = 50;  % samples; consecutive saccades closer than this are merged

        merged       = [];
        currentStart = blocks(1,1);
        currentEnd   = blocks(1,2);

        for i = 2:size(blocks, 1)
            nextStart = blocks(i,1);
            nextEnd   = blocks(i,2);

            if nextStart - currentEnd <= mergeGap
                % Close enough — extend the current saccade window
                currentEnd = max(currentEnd, nextEnd);
            else
                merged       = [merged; currentStart, currentEnd]; %#ok<AGROW>
                currentStart = nextStart;
                currentEnd   = nextEnd;
            end
        end

        saccades = [merged; currentStart, currentEnd];
    else
        saccades = blocks;
    end

end