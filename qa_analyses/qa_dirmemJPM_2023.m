% June 12, 2024. Updated code to be more than minimally useful. JPM
    % Changes: 
        % 1) Fixed sacc-aligned polar plot so that is actually shows
        % something now.
        % 2) Fixed black triangle plot so it now shows fitted 'best'
        % direction. WARNING: STILL NOT WORKING FOR SACC ALIGN
        % 3) changed color scheme from red-green to black-green.
        % 4) made linewidths larger in polar plots because I'm getting old
        % 5) added filename to plot


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

%%%%%%%%%%%%% MUST CHANGE THIS FOR YOUR PARAMS %%%%%%%%%%%%%
% The actual values here don't make a difference, but specify number of
% conditions in the order of appearance in xml file
% SETS{1} = saccade amplitudes
% SETS{2} = saccade angles
SETS = {[150,200], [0:45:315]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 3)
    sortcode = 1;
end

if (nargin < 4)
    h=figure;
    hb=figure;
else
    if (h > 0)
        f1 = figure(h);
        f1.Position = [100,50,1300,900];
        hb = h+1;
        f2 = figure(hb);
        f2.Position = [100,50,1300,900];
    end
end

a = qa_extract_data_1chan(filename,'mdir',sortcode,channel,SETS);

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

    for I=1:size(a.STIMEVENTS,1)
        if I == 1
            subplot(3,3,I+5)
        elseif I == 3 || I == 5
            subplot(3,3,I-1)
        elseif I == 4
            subplot(3,3,I-3)
        else % plotd = 2, 6-8
            subplot(3,3,I+1)
        end

        showPSTH(a.STIMEVENTS(I,:),xlimits,psthsmooth); box off;
      %  hold on
        overlayRaster(a.STIMEVENTS(I,:),xlimits,rastcol);
        xline(0, 'r-')
        set(get(gca,'Children'),'linewidth',psthline);

        % hard coded kludge; I==4
        if (I==4)
            xlabel('Aligned to Stim Onset. (s)');
            title ({filename; sprintf('sort code = %d',sortcode)}, 'Interpreter','none')
        end

        if (I<size(a.STIMEVENTS,1))

        else
            xlabel('Time (s)');
        end
    end
    set(gcf,'color','w')


    figure(hb)
    for I=1:size(a.STIMEVENTS,1)
        if I == 1
            subplot(3,3,I+5)
        elseif I == 3 || I == 5
            subplot(3,3,I-1)
        elseif I == 4
            subplot(3,3,I-3)
        else % plotd = 2, 6-8
            subplot(3,3,I+1)
        end
        
        overlayRaster(a.SACEVENTS(I,:),xlimits,rastcol);
        hold on
        showPSTH(a.SACEVENTS(I,:),xlimits,psthsmooth); box off;
        xline(0, 'r-')
        xline(0.26, 'r--')
        set(get(gca,'Children'),'linewidth',psthline);

        if (I==4)
            xlabel('Aligned to Sacc Onset. (s)');
            title ({filename; sprintf('sort code = %d',sortcode)}, 'Interpreter','none')
        end

        if (I<size(a.SACEVENTS,1))
        else
            xlabel('Time (s)');
        end
    end
    set(gcf,'color','w')
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


% Generate polar plots
if (h > 0)
   % figure(hb); clf;
   figure(h)
   subplot(3,3,5); 
    
 %   plim = max([nanmean(stimrate),nanmean(sacrate)]); 
    % can multiple * 1.1 to expand axes a bit
    
    % Plot stimulus onset aligned activity
 %   subplot(1,2,1); 
    rho = nanmean(stimrate);
    dst = sprintf('%0.2f',visds);
    dpt = sprintf('%0.2f',visdp);
    polarplot(deg2rad([theta 0]),[rho rho(1)],'ko-',...
        'markerfacecolor','k','linewidth',3)
    hold on
    polarplot(deg2rad([theta 0]),[rhoLst rhoLst(1)],'go--','LineWidth',2);
    polarplot(deg2rad([theta 0]),[rhoUst rhoUst(1)],'go--','LineWidth',2);   
    
    h1=polarplot(deg2rad(visdp),max(rho),'k^'); % Plot black triangle at best dir
    txd1 = get(h1,'ThetaData');
    tyd1 = get(h1,'RData');
    set(h1,'ThetaData',txd1(1),'RData',tyd1(1));
    set(h1,'markersize',10,'markerfacecolor','k'); hold off;
    title(['VisDir: ',dpt,', Sel: ',dst,', VMI: ',vmit]);
    prettyFig
    
    % Plot saccade-aligned responses
    figure(hb)
    subplot(3,3,5); 
    rho2 = nanmean(sacrate);
    dst = sprintf('%0.2f',sacds);
    dpt = sprintf('%0.2f',sacdp);
    polarplot(deg2rad([theta 0]),[rho2 rho2(1)],'ko-',...
        'markerfacecolor','k','linewidth',3)
    hold on
    polarplot(deg2rad([theta 0]),[rhoLsc rhoLsc(1)],'go--','LineWidth',2);
    polarplot(deg2rad([theta 0]),[rhoUsc rhoUsc(1)],'go--','LineWidth',2);

    h2=polarplot(sacdp,max(rho2),'k^'); % Plot black triangle at best dir
    txd2 = get(h2,'ThetaData');
    tyd2 = get(h2,'RData');
    set(h2,'ThetaData',txd2(1),'RData',tyd2(1));
    set(h2,'markersize',10,'markerfacecolor','k'); hold off;
    title(['SacDir: ',dpt,', Sel: ',dst,', STI: ',dsti]);
    prettyFig
    
    saveas(f1,sprintf('qa_figs/%s_unit%d_stimOnset.png',filename,sortcode),'png')
    saveas(f2,sprintf('qa_figs/%s_unit%d_saccOnset.png',filename,sortcode),'png')
end
end

function [all] = merge_raster(rasters, onemax)

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

