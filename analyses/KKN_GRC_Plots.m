%% GRC FIGURES

data_path = '/Users/kendranoneman/Data/dualhemi_unleashed';

%% individual session
sess = 'kendra_scrappy_0142a_g0';

load(fullfile(data_path,[sess '_unleashed.mat']), 'S')

%% histograms

imec0 = S.kilosort(1).clusters;
imec1 = S.kilosort(2).clusters;

imec0 = imec0(imec0.mdir1_Hz>1 & imec0.purs1_Hz>1,:);
imec1 = imec1(imec1.mdir1_Hz>1 & imec1.purs1_Hz>1,:);

imec0 = imec0(imec0.vis_sel_dir>0.1 | imec0.sac_sel_dir>0.1 | imec0.pur_sel_dir>0.1,:);
imec1 = imec1(imec1.vis_sel_dir>0.1 | imec1.sac_sel_dir>0.1 | imec1.pur_sel_dir>0.1,:);

f1a = figure;
f1a.Position = [100 100 1200 400];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
values = imec0.purs1_Hz;
histStyle(values,'left','ypos','number of clusters',[min(values)-10 max(values)+10],[0 300], 1, 10, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.purs1_Hz;
histStyle(values,'left','ypos','number of clusters',[min(values)-10 max(values)+10],[0 300], 1, 10, [0 0 0]./255, [100 100 100]./255)

%%

f1a = figure;
f1a.Position = [100 100 1200 400];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
values = imec0.vis_sel_dir;
histStyle(values,'left','ypos','number of clusters',[min(values)-1 max(values)+1],[0 300], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.vis_sel_dir;
histStyle(values,'left','ypos','number of clusters',[min(values)-1 max(values)+1],[0 300], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

%%

f1a = figure;
f1a.Position = [100 100 1200 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
values = imec0.VMI;
histStyle(values,'left','VMI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.VMI;
histStyle(values,'right','VMI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec0.SPI;
histStyle(values,'left','SPI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)

nexttile
values = imec1.SPI;
histStyle(values,'right','SPI','number of clusters',[-1 1],[0 100], 0.1, 0.8, [0 0 0]./255, [100 100 100]./255)


%% 

f2a = figure;
f2a.Position = [100 100 1200 800];
tl = tiledlayout(1,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
plot(imec0.VMI,imec0.SPI,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('SPI')
prettyFig;

nexttile
plot(imec1.VMI,imec1.SPI,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('SPI')
prettyFig;

%% 

f3a = figure;
f3a.Position = [100 100 1200 800];
tl = tiledlayout(2,2);
tl.TileSpacing = 'compact';
tl.Padding = 'tight';

nexttile
plot(imec0.SPI,imec0.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('SPI')
ylabel('ypos')
prettyFig;

nexttile
plot(imec1.SPI,imec1.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('SPI')
ylabel('ypos')
prettyFig;

nexttile
plot(imec0.VMI,imec0.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('ypos')
prettyFig;

nexttile
plot(imec1.VMI,imec1.y_pos,'ko')
hold on;
xline(0,'k--')
yline(0,'k--')
xlabel('VMI')
ylabel('ypos')
prettyFig;
