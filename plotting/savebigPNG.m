function savebigPNG(fignum, outputname)
set(gcf,'renderer','Painters') 

h=figure(fignum);set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,outputname,'-dpng','-r600')
end