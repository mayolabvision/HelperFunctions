function rgb = getColorFromValue(val, palette)
    if nargin < 2
        palette = '';
    end

    % Define endpoints
    switch upper(palette)
        case 'SPI'
            color2 = [242, 127, 57];  % SPI = -1 (pursuity, orange)
            color1 = [137, 110, 243];   % SPI = 1 (saccadey, purple)
        case 'VMI'
            color1 = [137, 110, 243];  % VMI = -1 (motor, purple)
            color2 = [10, 190, 155];  % VMI = 1 (visual, teal)
        case 'VPI'
            color1 = [242, 127, 57]; % VPI = -1 (pursuity, orange)
            color2 = [10, 190, 155]; % VPI = 1 (visual, teal)
        otherwise
            color1 = [255, 255, 255]; % white
            color2 = [0, 0, 0];       % black
    end

    val = max(min(val(:), 1), -1); % clip to [-1, 1]
    t = (val + 1) / 2;             % normalize to [0, 1]

    % Vectorized linear interpolation
    rgb = (1 - t) .* color1 + t .* color2;
    rgb = rgb / 255;
end