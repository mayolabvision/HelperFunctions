function prettyFig(varargin)
% PRETTYFIG Make figure look nice with thicker axes and ticks

ax = gca;  % current axes

% Ticks and font
set(ax, 'TickDir', 'out', ...
        'FontSize', 16, ...
        'LineWidth', 2);  % <-- thicker axes lines and tick marks

box off;

% Increase all text in figure
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'FontSize',14);

% Optionally set background color
if ~isempty(varargin)
    set(ax,'Color', varargin{1});
end

% Make figure background white
set(gcf,'Color','w');

% Make all plotted lines thicker (optional)
lines = findall(figureHandle,'type','line');
%set(lines,'LineWidth',2);

% Adjust figure size to match paper size
sz = get(gcf,'Position');
sz = sz(3:4); % width and height
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',sz);
end