function cmap = customColormap(color1, color2, N)
    % Inputs: color1, color2 as RGB triples in [0,255]; N = number of colors
    % Output: Nx3 colormap matrix scaled to [0,1]

    color1 = color1 / 255;
    color2 = color2 / 255;
    cmap = zeros(N, 3);
    for i = 1:3
        cmap(:, i) = linspace(color1(i), color2(i), N);
    end
end
