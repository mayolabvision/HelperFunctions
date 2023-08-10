% JPM

function savebigPDF(fignum, outputname)

% Test line added May 3 2019. Prevents large images from being switched
% from vectors to bitmap (e.g, lots of eye traces)
set(gcf,'renderer','Painters') 

h=figure(fignum);set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,outputname,'-dpdf','-r600')