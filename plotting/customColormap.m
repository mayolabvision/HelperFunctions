function cmap = customColormap(color1, color2, N)
    % Inputs: color1, color2 as RGB triples in [0,255]; N = number of colors
    % Output: Nx3 colormap matrix scaled to [0,1]

    if color2~=0
        color1 = color1 / 255;
        color2 = color2 / 255;
        cmap = zeros(N, 3);
        for i = 1:3
            cmap(:, i) = linspace(color1(i), color2(i), N);
        end
    else
        color1 = color1 / 255;
        % Create a brightness gradient from 0 (black) to 1 (white)
        t = linspace(0, 1, N)';
        
        % Convert midpoint to HSV to preserve hue/saturation
        mid_hsv = rgb2hsv(color1);
        
        % Define how brightness (value) varies: darker → mid → lighter
        % We'll go from 0.2*V to 1.2*V, clipped to [0,1]
        V = mid_hsv(3);
        V_scale = linspace(0.2*V, min(1.2*V,1), N)';
        
        % Construct new HSV values with constant hue/saturation but changing value
        custom_hsv = [repmat(mid_hsv(1), N, 1), ...  % constant hue
                       repmat(mid_hsv(2), N, 1), ...  % constant saturation
                       V_scale];                      % varying brightness
        
        % Convert back to RGB
        cmap = hsv2rgb(custom_hsv);
    end
end
