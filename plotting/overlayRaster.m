function overlayRaster(spikeTimes,timeRange,col)
%function overlayRaster(spikeTimes,timeRange,col)
%
% meant to be called right after showHist
%
%spikeTimes is a list of spikeTimes (in ms)
%  if spikeTimes is a cell array, all values inside will be pooled
%
% col is optional
%

if (nargin < 3)
    col = 'k';
end
  
numElements = size(spikeTimes,1);
numClasses = numel(spikeTimes)/numElements;

spikeTimes = reshape(spikeTimes,numel(spikeTimes),1);

if nargin > 1
  for i = 1:length(spikeTimes)
    spikeTimes{i} = spikeTimes{i}(spikeTimes{i} >= timeRange(1) & ...
                                  spikeTimes{i}<= timeRange(end));
  end
end

% necessary to scale the y-positions to overlay on showPSTH
yl = get(gca,'ylim');
yd = yl(2)-yl(1);
yi = yd/numClasses;
%yi = floor(yd/numClasses); % took this line out because of error when yi=0 (too many trials)
rs = yl(1):yi:yl(2); % raster steps
for i = 1:length(spikeTimes)
    h = line(repmat(spikeTimes{i}(:)',2,1),repmat([rs(i) rs(i+1)]',1,length(spikeTimes{i})));

    % old method from showRaster
    %h = line(repmat(spikeTimes{i}(:)',2,1),repmat([i-1 i]',1,length(spikeTimes{i})));    
    
    %if diff(timeRange) < 5
    %  set(h,'LineWidth',1)
    %end
    
    if ~isempty(h)
        set(h,'Color',col);%[mod(floor(i/122),2) 0 1-mod(floor(i/122),2)])
        set(h,'linewidth',0.5);
    end
    spikeTimes{i} = spikeTimes{i}(:);
end

axis tight


%for i = 1:numClasses-1
%  h = line(get(gca,'XLim'),[i*numElements i*numElements]);
%  set(h,'Color','k');
%  set(h,'LineWidth',2);
%end

%set(gca,'YTick',((1:numClasses)-.5)*numElements);
%set(gca,'YTickLabel',1:numClasses);

xlim(timeRange);
xlabel('Time (s)');