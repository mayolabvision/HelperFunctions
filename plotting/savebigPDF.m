% JPM w/ KKN additions

function savebigPDF(fignum, outputname, varargin)

p = inputParser;
addRequired(p, 'fignum');
addRequired(p, 'outputname', @(x) ischar(x) | isstring(x));
addOptional(p,'RESOLUTION', 600, @isnumeric) % DPI

p.parse(fignum, outputname, varargin{:});
fignum = p.Results.fignum;
outputname = p.Results.outputname;
RESOLUTION = p.Results.RESOLUTION;

h = figure(fignum); 
set(h,'Units','Inches');

% Test line added May 3 2019. Prevents large images from being switched
% from vectors to bitmap (e.g, lots of eye traces)
set(gcf,'renderer','Painters') 

pos = get(h,'Position');
set(h,'PaperPositionMode','Auto',...
      'PaperUnits','Inches',...
      'PaperSize',[pos(3), pos(4)])

resFlag = ['-r' num2str(RESOLUTION)];
print(h, outputname, '-dpdf', resFlag)

end