function rgb = getColorFromValue(val, palette)
    if nargin < 2
        palette = '';
    end

    % Define endpoints
    switch upper(palette)
        case 'SPI'
            color2 = [115, 30, 245];  % SPI = -1 (pursuity, purple)
            color1 = [245, 50, 30];   % SPI = 1 (saccadey, red)
        case 'VMI'
            color1 = [133, 255, 19];  % VMI = -1 (motor, green)
            color2 = [19, 133, 255];  % VMI = 1 (visual, blue)
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