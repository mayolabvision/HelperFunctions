clear
clc

addpath(genpath('/Users/kendranoneman/Packages')) % add nevUtils and HelperFunctions to path
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions')) 

DATAFOLDER = '/Users/kendranoneman/Data/SCRAPPY/raw';
DATAFILE = 'scrappy_0040_pursuit_behavOnly';

tbl = mayoFormat_trialDataTable(DATAFOLDER,DATAFILE);
correct_tbl = tbl(tbl.trialOutcome=='CORRECT',:);

%% Individual correct trials
speeds = sort(unique(tbl.pursuitSpeed));
angles = sort(unique(tbl.angle));

normalized_angles = angles / 360;
colors_rgb = hsv2rgb([normalized_angles, ones(length(angles), 2)]);
disp(colors_rgb);

for i = 1:height(correct_tbl)
    fig = figure;
    fig.Position = [100 100 1500 900];
    tl = tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'loose');

    x = 1:size(correct_tbl.eyePos{i},2);

    ax1 = nexttile;
    plot(x,correct_tbl.eyePos{i}(1,:),'linestyle','-','linewidth',2,'color',colors_rgb(correct_tbl.angle(i)==angles,:))
    hold on;
    plot(x,correct_tbl.eyePos{i}(2,:),'linestyle','--','linewidth',2,'color',colors_rgb(correct_tbl.angle(i)==angles,:))
    xlim([correct_tbl.TARG_ON(i)-100,correct_tbl.REWARD(i)])
    
    y_max = ylim; 

    xline(correct_tbl.TARG_ON(i), 'k', 'LineWidth', 1);
    text(correct_tbl.TARG_ON(i), y_max(2)+1, 'TARG ON', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    xline(correct_tbl.TARG_OFF(i), 'k', 'LineWidth', 1);
    text(correct_tbl.TARG_OFF(i), y_max(2)+1, 'TARG OFF', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

    xline(correct_tbl.REWARD(i), 'k', 'LineWidth', 1);
    text(correct_tbl.REWARD(i), y_max(2)+1, 'REWARD', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    ylabel('eye position [deg]')
    title(sprintf('%s',correct_tbl.trialName(i)),'interpreter','latex','fontsize',20,'fontweight','bold')
    subtitle(sprintf('speed = %d, angle = %d',correct_tbl.pursuitSpeed(i),correct_tbl.angle(i)))
    prettyFig;

    ax2 = nexttile;
    plot(x,correct_tbl.eyeVel{i}(1,:),'linestyle','-','linewidth',2,'color',colors_rgb(correct_tbl.angle(i)==angles,:))
    hold on;
    plot(x,correct_tbl.eyeVel{i}(2,:),'linestyle','--','linewidth',2,'color',colors_rgb(correct_tbl.angle(i)==angles,:))
    yline(correct_tbl.pursuitSpeed(i),'k--')
    text(correct_tbl.TARG_ON(i)-95, correct_tbl.pursuitSpeed(i)+1, sprintf('%d deg/s',correct_tbl.pursuitSpeed(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    xlim([correct_tbl.TARG_ON(i)-100,correct_tbl.REWARD(i)])
    %ylim([-(correct_tbl.pursuitSpeed(i)+15),(correct_tbl.pursuitSpeed(i)+15)])
    ylabel('eye velocity [deg/s]')
    prettyFig;

    ax3 = nexttile;
    plot(x,correct_tbl.eyeAcc{i}(1,:),'linestyle','-','linewidth',2,'color',colors_rgb(correct_tbl.angle(i)==angles,:))
    hold on;
    plot(x,correct_tbl.eyeAcc{i}(2,:),'linestyle','--','linewidth',2,'color',colors_rgb(correct_tbl.angle(i)==angles,:))
    xlim([correct_tbl.TARG_ON(i)-100,correct_tbl.REWARD(i)])
    ylabel('eye acceleration [deg/s^2]')
    prettyFig;

    linkaxes([ax1, ax2, ax3], 'x');
    xlabel(tl, 'Time [ms]','fontsize',18);
    
   

    pos = correct_tbl.eyePos{i};
    

end