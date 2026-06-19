% JPM w/ KKN additions, SVG version
function savebigSVG(fignum, outputname)
    p = inputParser;
    addRequired(p, 'fignum');
    addRequired(p, 'outputname', @(x) ischar(x) | isstring(x));
    p.parse(fignum, outputname);
    fignum = p.Results.fignum;
    outputname = p.Results.outputname;

    h = figure(fignum);
    set(h, 'Units', 'Inches');
    set(h, 'Renderer', 'Painters');   % keep it vector
    pos = get(h, 'Position');
    set(h, 'PaperPositionMode', 'Auto', ...
           'PaperUnits',        'Inches', ...
           'PaperSize',         [pos(3), pos(4)]);

    print(h, outputname, '-dsvg');
end