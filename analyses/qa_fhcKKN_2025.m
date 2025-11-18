function qa_fhcKKN_2025(filename,channel,sortcode)
% alignCode = code you want to align spikes to (x = 0)
% START_TRIAL = 1, STIM_ON = 10, TARG_ON = 70, SACCADE = 141, PURSUIT_TARG_ON = 31791
alignCode = 10; 

% condSplitBy = conditions you want to have as separate subplots (max 2)
% 'recColor', 'angle', 'distance', 'speed', 'jump'
condSplitBy = {'recColor'};

% x_limits = time window (in ms), aligned to alignCode, to plot spikes
% fr_window = time window (in ms) to plot shaded bar to highlight 
x_limits = [-100,200];
fr_window = [-50,50];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nev,~] = read_nev(filename);

trlStarts = find(nev(:,1)==0 & nev(:,2)==1);
trlStops = find(nev(:,1)==0 & nev(:,2)==255);

trldata = cell(numel(trlStops),2);
for t = 1:numel(trlStops)
   this_nev = nev(trlStarts(t):trlStops(t),:);
   
   trldig = this_nev(this_nev(:,1) == 0, :);
   if sum(trldig(:,2) == 150) == 0
       continue
   end
   trltext = char(trldig(trldig(:,2) >= 256 & trldig(:,2) < 512, 2) - 256)';
   
   trlspks = this_nev(this_nev(:,1) == channel & this_nev(:,2) == sortcode, 3);
   alignTime = trldig(trldig(:,2)==alignCode,3);
   trlspks = (trlspks - alignTime(1))*1000;
   
   trldata{t,1} = trltext;
   trldata{t,2} = trlspks;
end

T = cell2table(trldata(sum(cellfun(@(q) isempty(q), trldata, 'uni', 1),2)==0,:), 'VariableNames', {'conditions','spiketimes'});

if isempty(T)
    fprintf('Error: No data found for that sortcode\n')
    return 
end

conditions_split = cellfun(@(x) strsplit(x, ';'), T.conditions, 'uni', 0);
first_condition = conditions_split{1}; 
keys = cellfun(@(x) strsplit(x, '='), first_condition, 'uni', 0);
condition_names = cellfun(@(x) x{1}, keys, 'uni', 0);
condition_names = condition_names(1:end-1);

condition_values = cell2mat(cellfun(@(v) cellfun(@(q) str2double(regexp(q, '\d+', 'match')), v(1:end-1), 'uni', 1), conditions_split, 'uni', 0));

for c = 1:numel(condition_names)
   T.(condition_names{c}) = condition_values(:,c);
end

fig = figure;
fig.Position = [100 100 1100 900];

if isempty(condSplitBy)
    raster_sdf(T.spiketimes, 'FR_WINDOW', fr_window, 'TIME_WINDOW', x_limits)
else
    conds = cellfun(@(q) sort(unique(T.(q))), condSplitBy, 'uni', 0);

    plt_num = 0;
    y_lims = [];
    if numel(conds) == 1
         for r = 1:numel(conds{1})
             plt_num = plt_num + 1;
             subplot(numel(conds{1}), 1, plt_num)
             these_trls = T(T.(condSplitBy{1})==conds{1}(r),:);

             raster_sdf(these_trls.spiketimes, 'FR_WINDOW', fr_window, 'TIME_WINDOW', x_limits)

             yyaxis left;
             ax = gca;
             y_lims = [y_lims; ax.YLim];

             title(sprintf('%s = %d',string(condSplitBy{1}),conds{1}(r)),'FontSize',14)
         end
    else
        for r = 1:numel(conds{1})
            for c = 1:numel(conds{2})
                plt_num = plt_num + 1;
                subplot(numel(conds{1}), numel(conds{2}), plt_num)
                these_trls = T(T.(condSplitBy{1})==conds{1}(r) & T.(condSplitBy{2})==conds{2}(c),:);

                raster_sdf(these_trls.spiketimes, 'FR_WINDOW', fr_window, 'TIME_WINDOW', x_limits)

                yyaxis left;
                ax = gca;
                y_lims = [y_lims; ax.YLim];

                title(sprintf('%s = %d, %s = %d',string(condSplitBy{1}),conds{1}(r),string(condSplitBy{2}),conds{2}(c)))  
            end

        end
    end

    global_y_lim = [min(y_lims(:,1)), max(y_lims(:,2))];

     for p = 1:plt_num
        if numel(conds) == 1
            subplot(numel(conds{1}), 1, p)
            yyaxis left;
            ylim(global_y_lim);
        else
            subplot(numel(conds{1}), numel(conds{2}), p)
            yyaxis left;
            ylim(global_y_lim);
        end
     end
end

end

