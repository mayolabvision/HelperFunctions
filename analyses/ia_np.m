%% 
load('/Volumes/home/DATA/kendra_scrappy_0127a_g0/kendra_scrappy_0127a.mat')

%% b

[frs,bin_edges,xvals,yvals] = format_tableToRFMap(S.rfmp1.data,FIRST_BIN,BIN_WIDTH,BIN_STEP,NBINS);

%% 
fig_path = '/Volumes/home/DATA/kendra_scrappy_0127a_g0/figs/rfmp/unit_heatmaps';
filename = 'kendra_scrappy_0127a';

for unit=14:length(frs)
    f2a = figure('Visible','off');
    f2a.Position = [100 100 1800 900];
    tl = heatMap_rfOverTime(frs{unit},'BIN_EDGES',bin_edges, 'INTERP', false,'X_VALS',pix2deg(xvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)), 'Y_VALS',pix2deg(yvals,S.rfmp1.params.screenDistance(1),S.rfmp1.params.pixPerCM(1)));
    
    title(tl,sprintf('%s --- unit %d',filename,unit),'fontsize',20,'interpreter','none')
    %subtitle(tl,sprintf('%s (ripChan = %d, depth = %2.3f mm)',chan_name, S.channels.ripChan_num(good_chans(unit)), chan_depth),'fontsize',16,'interpreter','none')
    
    print(f2a, fullfile(fig_path, sprintf('unit%3d.png', unit)), '-dpng', '-r200');
    fprintf(sprintf('\n----Unit %.2d complete----',unit))

    blah = 1;
end

fprintf('\n----------------------\n')