function rgb = radial_colormap(dir, varargin)
% radial_colormap returns an RGB color (0–1) for a given direction (0–359 deg)
% Usage:
%   rgb = radial_colormap(dir)
%   rgb = radial_colormap(dir, 'SHIFT_DEG', 90)
%
% dir can be scalar or vector. Output is Nx3.

    % ------------------------------
    % Parse inputs
    % ------------------------------
    p = inputParser;
    addRequired(p, 'dir', @isnumeric);
    addParameter(p, 'SHIFT_DEG', 0, @isnumeric);
    parse(p, dir, varargin{:});
    dir = p.Results.dir;
    shift_deg = p.Results.SHIFT_DEG;

    % ------------------------------
    % Define canonical 8-color wheel
    % ------------------------------
    key_angles = [0 45 90 135 180 225 270 315];

    key_colors = [
        1.00 0.00 0.00   % 0°   red
        1.00 0.55 0.00   % 45°  orange
        1.00 1.00 0.00   % 90°  yellow
        0.00 0.80 0.00   % 135° green
        0.20 0.80 1.00   % 180° light blue
        0.00 0.00 0.90   % 225° dark blue
        0.55 0.00 0.80   % 270° purple
        1.00 0.40 0.70   % 315° pink
    ];

    % ------------------------------
    % Rotate and sort the key angles
    % ------------------------------
    shifted_angles = mod(key_angles + shift_deg, 360);
    [shifted_angles, idx] = sort(shifted_angles);
    key_colors = key_colors(idx, :);

    % ------------------------------
    % Make extended (circular) lookup
    % ------------------------------
    angles_ext = [shifted_angles, shifted_angles + 360];      % 1 x 16
    colors_ext = [key_colors; key_colors];                    % 16 x 3

    % ------------------------------
    % Prepare query angles
    % ------------------------------
    dir_mod = mod(dir, 360);          % map into [0,360)
    dirq = dir_mod(:);                % column vector queries

    % If the query is less than the first shifted angle, shift it up by 360
    % so it falls within the extents of angles_ext (vectorized)
    minAng = angles_ext(1);
    mask = dirq < minAng;
    dirq(mask) = dirq(mask) + 360;

    % ------------------------------
    % Interpolate R, G, B separately (no extrapolation needed)
    % ------------------------------
    rgb = zeros(numel(dirq), 3);
    for c = 1:3
        rgb(:, c) = interp1(angles_ext, colors_ext(:, c), dirq, 'linear');
    end

    % ------------------------------
    % Safety: clip into [0,1]
    % ------------------------------
    rgb = max(0, min(1, rgb));

    % Restore shape: if input was row, return row(s) accordingly
    if isrow(dir)
        rgb = reshape(rgb, size(dir,2), 3);   % return Nx3 where N = numel(dir)
    end
end