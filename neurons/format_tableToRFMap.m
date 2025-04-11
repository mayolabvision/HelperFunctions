function [all_FRs,bin_edges,xvals,yvals] = format_tableToRFMap(T,first_bin,bin_width,bin_step,nbins)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    stimuli = cell2mat(vertcat(T.conditions{:}));
    xvals = sort(unique(stimuli(:,1)));
    yvals = sort(unique(stimuli(:,2)));

    bin_edges = arrayfun(@(x) [(first_bin + ((x*bin_step)-bin_step)),(first_bin + ((x*bin_step)-bin_step) + bin_width)], 1:nbins, 'UniformOutput', false);

    num_channels = length(T.spiketimes_imec0{1});

    all_FRs = cell(num_channels,1);
    for unit = 1:num_channels
        FRs = cell(length(yvals),length(xvals),length(bin_edges));
        for bin = 1:length(bin_edges)
            for trial = 1:height(T)
                stim_ons = T.STIM_ON{trial};
                spks = T.spiketimes_imec0{trial}{unit};
                %nets = T.net_labels{trial,unit};
                for stim = 1:length(stim_ons)
                    aligned_spks = spks-stim_ons(stim);
                    %spk_hz = sum(aligned_spks>bin_edges{bin}(1) & aligned_spks<=bin_edges{bin}(2) & nets>GAMMA)*(1000/range(bin_edges{bin}));
                    spk_hz = sum(aligned_spks>bin_edges{bin}(1) & aligned_spks<=bin_edges{bin}(2))*(1000/range(bin_edges{bin}));
                    
                    if isempty(FRs{T.conditions{trial}{stim}(2)==yvals,T.conditions{trial}{stim}(1)==xvals,bin})
                        FRs{T.conditions{trial}{stim}(2)==yvals,T.conditions{trial}{stim}(1)==xvals,bin} = spk_hz;
                    else
                        FRs{T.conditions{trial}{stim}(2)==yvals,T.conditions{trial}{stim}(1)==xvals,bin} = [FRs{T.conditions{trial}{stim}(2)==yvals,T.conditions{trial}{stim}(1)==xvals,bin}, spk_hz];
                    end
                end
            end
        end
        all_FRs{unit} = FRs;
        fprintf(sprintf('\n----Unit %.2d complete----',unit))
    end
    fprintf('\n----------------------\n')
end