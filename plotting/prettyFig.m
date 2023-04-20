function prettyFig(varargin)

set(gca, 'TickDir', 'out','FontSize',16)
box off
%axis square
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',16) % enlarge all text in figure to size 14
%set(findall(figureHandle,'type','line'),'LineWidth',2) % make all lines thicker to width 2
if length(varargin)~=0
    set(gca,'color',varargin{1})
end
set(gcf,'color','w')


sz = get(gcf,'Position');
sz = sz(3:4); % the last two elements are width and height of the figure
set(gcf,'PaperUnit','points'); % unit for the property PaperSize
set(gcf,'PaperSize',sz);