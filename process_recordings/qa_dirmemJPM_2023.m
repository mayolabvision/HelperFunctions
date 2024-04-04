% May 5, 2023. My triumphant return to and updating of this code that 
% I apparently wrote over 10 years ago. Goal is to make the plots pretty
% and immediately legible, unlike whatever pile of shit that is currently
% being used in the lab.  Apparently, Matlab changed polplot --> polar -->
% polarplot

function [visdp, visds, sacdp, sacds, vmi, STI, baseFR] = qa_dirmemJPM_2023(filename,channel,sortcode,h)
% function [visdp visds sacdp sacds vmi] = qa_dirmem(filename,sortcode,h)
%
% qa_dirmemJPM('Ya230504_s74n2_dirmem_0004.nev', 257, 1,8);

% sortcode is optional list of valid sort codes (default is 1)
% 
% h is an optional figure handle. If h is -1, this doesn't plot and
% instead just spits out the results

if (nargin < 3)
    sortcode = 1;
end

if (nargin < 4)
    h=figure;
    hb=figure;
else
    if (h > 0)
        figure(h);
        hb = h+1;
        figure(hb);
    end
end

a = qa_extract_data_1chan(filename,sortcode,channel); 

%%%%%%%
% Generate two raster/PSTH plots - stim-aligned and sac-aligned
%%%%%%%
xlimits = [-.3 .5];
psthsmooth = 20;
rastcol = [0.7 0.7 0.7];
psthline = 1.5;

sigma = 10;  % For spike density functions.
%

if (h > 0)
    figure(h); clf;
    % stim-aligned first
    pcount = [1:2:length(a.CND)*2];
    
    for I=1:size(a.STIMEVENTS,1)
        subplot(size(a.STIMEVENTS,1),2,pcount(I));
        showPSTH(a.STIMEVENTS(I,:),xlimits,psthsmooth); box off;
        rasters = showPSTH(a.STIMEVENTS(I,:),xlimits,psthsmooth); % JPM added
        set(get(gca,'Children'),'linewidth',psthline);
        if I==1  % if-statement added by JPM, June 25 2013
            ylabel(['CND: ',num2str(I), ' RIGHT']);
        elseif I==3
            ylabel(['CND: ',num2str(I), ' UP']);
        else
            ylabel(['CND: ',num2str(I)]);
        end
        if (I==1)
            title('Aligned to Stim Onset');
        end
        if (I<size(a.STIMEVENTS,1))
            set(gca,'xticklabel','');
        else
            xlabel('Time (s)');
        end
    end

    % sac-aligned
    pcount = [2:2:length(a.CND)*2];
    for I=1:size(a.SACEVENTS,1)
        subplot(size(a.SACEVENTS,1),2,pcount(I));
        showPSTH(a.SACEVENTS(I,:),xlimits,psthsmooth); box off;  
        set(get(gca,'Children'),'linewidth',psthline);    
        %ylabel(['CND: ',num2str(I)]);
        if (I==1)
            title('Aligned to Saccade Onset');
        end
        if (I<size(a.SACEVENTS,1))
            %set(gca,'xticklabel','');
        else
            xlabel('Time (s)');
        end
    end
    
    % set same ylimits on both axes
    pc = get(gcf,'children');
    for I=1:length(pc)
        yl(I,:) = get(pc(I),'ylim');
    end
    for I=1:length(pc)
        set(pc(I),'ylim',[0 max(yl(:,2))]);        
    end
    
    % overlay rasters for stim-aligned and sac-aligned
    pcount = [1:2:length(a.CND)*2]; 
        holdrasters=[];
    for I=1:size(a.STIMEVENTS,1)
        subplot(size(a.STIMEVENTS,1),2,pcount(I));
        ylim([min(yl(:,1)) max(yl(:,2))]);
        hold on
        overlayRaster(a.STIMEVENTS(I,:),xlimits,rastcol);
        plot([0 0],[0 max(yl(:,2))],'r');
        hold off
    end
    
    pcount = [2:2:length(a.CND)*2];    
    for I=1:size(a.SACEVENTS,1)
        subplot(size(a.SACEVENTS,1),2,pcount(I));
        ylim([min(yl(:,1)) max(yl(:,2))]);    
        hold on
        overlayRaster(a.SACEVENTS(I,:),xlimits,rastcol);
        plot([0 0],[0 max(yl(:,2))],'r');
        hold off
    end
end


%%%% calculate some statistics
% spike counting windows
stimwin = [0.05 0.15];
sacwin = [-0.05 0.05];

stimrate = nan(size(a.STIMEVENTS,1),size(a.STIMEVENTS,2));
for I=1:size(a.STIMEVENTS,1)
    tdat = squeeze(a.STIMEVENTS(I,1:a.RPTS(I)));
    for J=1:size(a.STIMEVENTS,2)
        if J > a.RPTS(I)
            break
        end
        stimrate(I,J) = length(find(tdat{J}>=stimwin(1) & tdat{J}<stimwin(2)));
    end

end
stimrate = stimrate' ./ diff(stimwin);

sacrate = nan(size(a.SACEVENTS,1),size(a.SACEVENTS,2));
for I=1:size(a.STIMEVENTS,1)
    tdat = squeeze(a.SACEVENTS(I,1:a.RPTS(I)));
    for J=1:size(a.SACEVENTS,2)
        if J > a.RPTS(I)
            break
        end
        sacrate(I,J) = length(find(tdat{J}>=sacwin(1) & tdat{J}<sacwin(2)));
    end
end
sacrate = sacrate' ./ diff(sacwin);

theta = 0:360/length(a.CND):360; theta(end)=[];

% JPM Aug 15, 2012 modifying previous version of qa_dirrem to
% permute/bootstrap orig firing rates 

% Generate randomized index of stimrate values, WITH REPLACEMENT
shuffles = 1000;
rhoPst=[]; rhoPsc=[];

for sh=1:shuffles
    randind=randi( (size(stimrate,1)*size(stimrate,2)), size(stimrate,1), size(stimrate,2) );
    permutedStimrate = stimrate(randind);
    rhoPst = [rhoPst; nanmean(permutedStimrate)];
    randind=randi( (size(sacrate,1)*size(sacrate,2)), size(sacrate,1), size(sacrate,2) );
    permutedSacrate = sacrate(randind);
    rhoPsc = [rhoPsc; nanmean(permutedSacrate)];
end

sorted_rhoPst=sort(rhoPst);
rhoLst = sorted_rhoPst(shuffles*.05,:); % 95% lower confidence interval
rhoUst = sorted_rhoPst(shuffles-(shuffles*.05),:); % 95% upper confidence interval

sorted_rhoPsc=sort(rhoPsc);
rhoLsc = sorted_rhoPsc(shuffles*.05,:); % 95% lower confidence interval
rhoUsc = sorted_rhoPsc(shuffles-(shuffles*.05),:); % 95% upper confidence interval


% calculate visuomotor index, where 1 is visual and -1 is motor
visual = nanmean(nanmean(stimrate));
motor = nanmean(nanmean(sacrate));
vmi = (visual - motor) / (visual + motor);
vmit = sprintf('%0.2f',vmi);


% BEGIN STI (sustained-transient index) code, from WhiteMunoz09
% find transient peak in visual response 50-150 ms after stimOn (MayoSommer 2013)
indexx=cellfun(@(x) x>=min(xlimits) & x<max(xlimits), a.STIMEVENTS, 'uni', 0); % index of spike times w/i relevant limits
allspktimes = cellfun(@(x, y) x(y), a.STIMEVENTS,indexx, 'uni', 0); % keep spike times that occur w/i relevant limits
spks2=reshape(allspktimes,numel(allspktimes),1);
spks2=cellfun(@(x) round((x*10000)+3000), spks2, 'uni', 0); % HARD CODED VALUES**

holdrasters = zeros(numel(allspktimes), range(xlimits*10000) ); % eg 80 (10 trials x 8 directions) x 8000 ms

for i=1:size(holdrasters,1)
    holdrasters(i,spks2{i,1}(spks2{i,1}>0))=1;
end

sdfstim = spike_density( merge_raster( holdrasters ), sigma ); % StimOn at 3000

peakfr = max(sdfstim(3500:4500));
peaklatency = find(sdfstim(3500:4500)==peakfr) + 499; % in ms relative to stimOn
halfmaxlatency = find(sdfstim(3000:(3000+peaklatency))<(peakfr/2), 1, 'last'); % time to half-max relative to stimOn, in 

sumrast = sum(holdrasters,1);

spksBASE = sum(sumrast(2000:3000)); % baseline total # spikes= spikes -100 to 0 ms before StimOn
baseFR = (spksBASE / numel(allspktimes) ) * 10; % convert from spikes per 100 ms to spikes per 1000 ms

% remember, stimOn occurs at 3000 in sdfstim and sumrast
spksTRNS = sum(sumrast(3000+halfmaxlatency:3500+halfmaxlatency)); % visual transient= # spikes 0-50 ms after time to half-max
spksSUST = sum(sumrast(3500+halfmaxlatency:4000+halfmaxlatency)); % sustained= # spikes 50-100 ms after time to half-max
STI = (spksTRNS-spksSUST)/(spksTRNS+spksSUST); % Positive values indicate transient larger than sustained
dsti = sprintf('%0.2f',STI);
% END STI



% calculate tuning preferences
%theta = 0:360/length(a.CND):360; theta(end)=[];
[visds, visdp] = tuningbias(theta,nanmean(stimrate));
[sacds, sacdp] = tuningbias(theta,nanmean(sacrate));

if (h > 0)
    % polar plots
    figure(hb); clf;
    
 %   plim = max([nanmean(stimrate),nanmean(sacrate)]); 
    % can multiple * 1.1 to expand axes a bit
    
    subplot(1,2,1);
    rho = nanmean(stimrate);
    dst = sprintf('%0.2f',visds);
    dpt = sprintf('%0.2f',visdp);
    polarplot(deg2rad([theta 0]),[rho rho(1)],'ro-',...
        'markerfacecolor','r','linewidth',2)
    hold on
    polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'go--');
    polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'go--');   
    %set(hP,'markerfacecolor','g');  % JPM
    h2=polarplot(deg2rad(visdp),max(rho),'k^');
%     txd = get(h2,'xdata');
%     tyd = get(h2,'ydata');
%     set(h2,'xdata',txd(1),'ydata',tyd(1));
%     set(h2,'markersize',10,'markerfacecolor','k'); hold off;
    title(['VisDir: ',dpt,', Sel: ',dst,', VMI: ',vmit]);
    prettyFigSGL
    
    subplot(1,2,2);
    rho = nanmean(sacrate);
    dst = sprintf('%0.2f',sacds);
    dpt = sprintf('%0.2f',sacdp);
    h1=polplot(theta,rho,'ro-'); hold on;
    set(h1,'markerfacecolor','r');
    polplot(theta,rhoLsc,'go--'); hold on; % JPM
    polplot(theta,rhoUsc,'go--'); hold on; % JPM   
    h2=polplot(sacdp,max(rho),'k^');
    txd = get(h2,'xdata');
    tyd = get(h2,'ydata');
    set(h2,'xdata',txd(1),'ydata',tyd(1));
    set(h2,'markersize',10,'markerfacecolor','k'); hold off;
    title(['SacDir: ',dpt,', Sel: ',dst,', STI: ',dsti]);
    prettyFigSGL
end
end

function [all] = merge_raster( rasters, onemax)

if nargin < 2
    onemax = 0;
end

sz = size( rasters );
if sz(1) <2
    all = rasters;
else
    all = sum(rasters);
end

if onemax
    f = find(all > 1);
    all(f) = 1;
end
end

