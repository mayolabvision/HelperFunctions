% Dec 20, 2024. added some cursory runNASnet functionality with the help of
% Deepa


function qa_rfmapMayoLab (nevname, savepdf, runnasnet)

% % qa = 'quick analysis'

% qa code for Mayo Lab (that Kendra will inevitably discard or entirely
% rewrite) for plotting RF maps using rfmap task. Currently discared the
% 'extra' trials where fewer stimuli are shown. But those few trials should
% be added back later. 

% for now, assuming equally spaced grid of 16 x 12 stimuli centered in
% screen

% TRIAL CODES
% 1 = START_TRIAL
% 2 = FIX_ON
% 140 = FIXATE (start of fixation)
% 10 = STIM_ON
% 40 = STIM_OFF
% 4 = FIX_MOVE
% 150 = CORRECT
% 3 = FIX_OFF
% 5 = REWARD

% DEPENDENT FUNCTIONS
% readNEV
% nev2dat
% savebigPDF

if nargin < 3
    runnasnet = 0;
    if nargin < 2
        savepdf = 0;
    end
end

if runnasnet == 1
    [slabel,sorted_spikes,net_labels,waveforms] = ...
        runNASNet(['/Users/jpm/Desktop/Dropbox/Dropbox/MAYO_LAB/data_files/Scrappy/' nevname],0.2, ...
        'netFolder', '/Users/jpm/Documents/GitHub/spikesort/nasnet/networks','netname', 'UberNet_N50_L1');
else
    disp 'USING ALL SPIKES ON CHANNEL, NO SORTING'
    nev = readNEV(nevname);
end



% extract the trial codes
dat2 = nev2dat(nevname,'readNS5',true,'convertEyes',true);

% matrix of image onset times, trials (rows) x images per trial (cols)
% iot3 = num2cell ( cellfun(@(x,y) x+(trlleng*(y-1)), imOntimes2, repmat(rowind, 1, 24))) ; 

trials = size(dat2, 2);

tt = struct2table(dat2);

stimOn_times_TrialsByImages = [];
xvals = [];
yvals = [];

for trialnum = 1:trials
    % trialnum % uncomment for debugging purposes

    % logical test if trialcodes for this trial contains 150 (equals 1)
    if sum(tt.trialcodes{trialnum,1}(:,2)==150)==1

        % logical index of stimulus onset times (code 10) in the trial
        ind_ons = tt.trialcodes{trialnum,1}(:,2)== 10;


        if length(tt.trialcodes{trialnum,1}(ind_ons,3)) == 11 % 10 or 11?

            % stimulus onset times in the trial, in Ripple time values. should
            % always be 11 values (for now)
            stimOn_times_TrialsByImages = [stimOn_times_TrialsByImages tt.trialcodes{trialnum,1}(ind_ons,3)];

            sup=dat2(trialnum).text; % store long disgusting list of characters that are x and y stim positions

            indxposstart = strfind(sup, 'xpos=')+5;
            indxposend = strfind(sup, ';ypos')-1;

            indYposstart = strfind(sup, 'ypos=')+5;
            indYposend = [ strfind(sup, ';xpos')-1 length(sup)];


            for jj = 1:length(indxposend)


                if isnan ( str2double(sup(indYposstart(jj):indYposend(jj))))
                    tempxvals(jj) = str2double(sup(indxposstart(jj):indxposend(jj))) ;
                    tempyvals(jj) = str2double(sup(indYposstart(jj):indYposend(jj)-1)) ;
                    
                else
                    tempxvals(jj) = str2double(sup(indxposstart(jj):indxposend(jj))) ;
                    tempyvals(jj) = str2double(sup(indYposstart(jj):indYposend(jj))) ;

                end
            end

            xvals = [ xvals; tempxvals];
            yvals = [ yvals; tempyvals];
            clear tempxvals tempyvals tempvalsBOTH

        end

    end
end



possible_stim_XYlocations = table2array ( combinations ( -495:90:495, -495:90:495 ) );


fig = figure;
fig.Units = 'inches';
fig.Position = [ 6 10 18 12];

% % PLOTTING BELOW
for probe_chan = 1:32 % size(allfns,1) % go thru each unit, and make a plot for each unit

    if runnasnet == 0
        allspks_thischan = nev(nev(:,1)==probe_chan, 3);
    else
        allspks_thischan = sorted_spikes(sorted_spikes(:,1)==probe_chan & sorted_spikes(:,2)==1, 3);
    end

    CLIMS=[];
    %  for w = 10:10:240
    for lat = 0.010:0.010:0.240  % go thru multiple spike windows.  in UNITS??

        frperloc = [];

        for dur = 0.050 % in UNITS?

         %   spks_ImageDims = repmat({allspks_thischan}, size(stimOn_times_TrialsByImages));

            stimOn_times_TrialsByImagesCELL = num2cell(stimOn_times_TrialsByImages);

            spikewin = [lat lat+dur];

            allspks_for_all_images = repmat({allspks_thischan}, size(stimOn_times_TrialsByImages)); % put all of unit's spikes in single cells to

            % this yields the summed number of spikes (per spikewin) in a matrix of Trial x Image
            sumspkperim = cellfun(@(x, y) sum(x > (y+spikewin(1)) & x <= (y+spikewin(2)) ), allspks_for_all_images, stimOn_times_TrialsByImagesCELL, 'uni', 0);
            sumspkperim = transpose (sumspkperim);


            for s = 1:144 % go thru each possible stimulus location

                ind_stimLocOnTrials = xvals==( possible_stim_XYlocations(s,1) ) & yvals== possible_stim_XYlocations(s,2);

                % safety check
                %if 

                frperloc(s) = mean (cell2mat(sumspkperim(ind_stimLocOnTrials)));

            end

            subplot(4,6,round(lat*100));
            imagesc(reshape(frperloc, 12, 12));
            colormap gray
            hcb = colorbar;

            title (['ch', num2str(probe_chan)  ' ' num2str(spikewin(1)), ':', num2str(spikewin(2)) ])
            prettyFig
            CLIMS = [CLIMS get(gca, 'clim')];
        end

    end


    % scale subplots
    for sp = 1:24
        subplot(4,6,sp)
        set( gca, 'clim', [min(CLIMS) max(CLIMS)] )
    end
 %   disp (['This is ', allfns{u}])
  %    pause

   if savepdf == 1
       if runnasnet == 0
       savebigPDF(1, [nevname(1:end-4), '_ch', num2str(probe_chan), '.pdf'] )
       elseif runnasnet == 1
        savebigPDF(1, [nevname(1:end-4), '_ch', num2str(probe_chan), '_SORTED.pdf'] )
       end
   else
       pause
   end

    clf


end % end of loop thru spike channels