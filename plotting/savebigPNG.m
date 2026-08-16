function savebigPNG(fignum, outputname, resolution)

% --- Default resolution ---
if nargin < 3 || isempty(resolution)
    resolution = 300;
end

% --- Get figure handle safely ---
h = figure(fignum);

% --- Use Painters renderer for vector-friendly export ---
set(h,'Renderer','painters');
set(h,'Units','Inches');

% --- Match paper size to on-screen size ---
pos = get(h,'Position');
set(h,'PaperPositionMode','auto', ...
      'PaperUnits','Inches', ...
      'PaperSize',[pos(3), pos(4)]);

% --- Print at requested resolution ---
print(h, outputname, '-dpng', ['-r' num2str(resolution)]);

% --- Close figure safely (only if still valid) ---
if isgraphics(h)
    close(h);
end

end