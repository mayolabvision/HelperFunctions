function rgb_color = radial_colormap(dir, varargin)
    % radial_colormap returns RGB color(s) for a given direction(s) in degrees
    %
    % Usage:
    %   rgb = radial_colormap(dir)                   % default shift: 0
    %   rgb = radial_colormap(dir, 'SHIFT_DEG', 90)
    %
    % Inputs:
    %   dir - direction(s) in degrees (scalar or vector)
    %   'SHIFT_DEG' - optional name-value pair (default: 0)
    %
    % Output:
    %   rgb_color - Nx3 matrix of RGB values in [0, 1] for each input direction

    p = inputParser;
    addRequired(p, 'dir', @isnumeric);
    addParameter(p, 'SHIFT_DEG', 0, @isnumeric);

    parse(p, dir, varargin{:});
    dir = p.Results.dir;
    shift_degrees = p.Results.SHIFT_DEG;

    % Define the key colors as RGB values (normalized to the range [0, 1])
    key_colors = [
        251, 176, 58;
        240, 90, 37;
        236, 29, 35;
        254, 1, 0;
        210, 20, 90;
        147, 39, 142;
        102, 46, 145;
        46, 49, 145;
        5, 168, 157;
        0, 146, 68;
        217, 224, 34;
        253, 237, 34
    ] / 255;

    % Corresponding positions (in degrees) of the key colors
    key_positions = [1, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330];

    % Rotate positions
    key_positions = mod(key_positions - shift_degrees, 360);

    % Sort positions and corresponding colors
    %[key_positions, sort_idx] = sort(key_positions);
    %key_colors = key_colors(sort_idx, :);

    % Preallocate the colors cell array with 359 cells
    colors = cell(1, 359);

    % Loop through each of the 359 positions
    for i = 1:359
        % Find the two nearest key colors to interpolate between
        idx1 = find(key_positions <= i, 1, 'last');
        idx2 = find(key_positions >= i, 1, 'first');
        
        % Handle edge cases for circular interpolation
        if isempty(idx1)
            idx1 = length(key_positions);
        end
        if isempty(idx2)
            idx2 = 1;
        end
        
        % Get the positions and colors for interpolation
        pos1 = key_positions(idx1);
        pos2 = key_positions(idx2);
        color1 = key_colors(idx1, :);
        color2 = key_colors(idx2, :);
        
        % Perform linear interpolation between the two colors
        if pos1 ~= pos2
            t = (i - pos1) / (pos2 - pos1);
            interpolated_color = (1 - t) * color1 + t * color2;
        else
            interpolated_color = color1;  % No interpolation needed
        end
        
        % Store the interpolated color in the cell array
        colors{i} = round(interpolated_color * 255);  % Convert back to [0, 255] range
    end

    % Output the RGB color for the given 'dir' (ensure dir is in valid range)
    dir = mod(dir-1, 359) + 1;  % Ensure 1 <= dir <= 359 for circular indexing
    rgb_color = colors{dir}./255;    % Return the RGB color at the index 'dir'
end
