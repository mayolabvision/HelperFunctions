function [all_FRs, bin_edges, xvals, yvals] = format_tableToRFMap(tbl, varargin)
    % Computes peristimulus time histogram (PSTH) binned firing rates for each stimulus condition in a 
    % table and returns data structured for making RF mapping plots
    %
    % INPUTS:
    %   tbl           - A table where each row is a trial. Required columns:
    %                   - 'conditions': each cell contains an Nx2 matrix of stimulus positions [x, y]
    %                   - 'spiketimes_imecX': where X is the IMEC probe number (e.g., spiketimes_imec0),
    %                     each cell is a 1xN cell array of spike times per unit (ms)
    %                   - 'STIM_ON': each cell is a vector of stimulus onset times for that trial, 
    %                      where time = 0 is trial onset (ms)
    %
    % Optional Name-Value Pair Arguments:
    %   'IMEC'        - Integer specifying which probe's spike times to use (default = 0)
    %   'FIRST_BIN'   - Time (ms) of the first bin relative to stimulus onset (default = 0)
    %   'BIN_WIDTH'   - Width (ms) of each time bin (default = 50)
    %   'BIN_STEP'    - Step size (ms) between successive bins (default = 10)
    %   'N_BINS'      - Number of bins to compute (default = 10)
    %
    % OUTPUTS:
    %   all_FRs   - Cell array where each cell contains a 3D cell matrix of firing rates:
    %               rows = unique y-values of stimulus, cols = unique x-values,
    %               depth = time bins
    %   bin_edges - Cell array of time bin edges used to compute firing rates
    %   xvals     - Sorted list of unique x stimulus positions
    %   yvals     - Sorted list of unique y stimulus positions
    %
    % EXAMPLE USAGE:
    %   [all_FRs, bin_edges, xvals, yvals] = format_tableToRFMap(tbl, 'IMEC', 1, 'BIN_WIDTH', 100);
    %
    %   This would extract spike times from the column 'spiketimes_imec1', compute firing rates
    %   in 100ms bins, and return a structure suitable for RF visualization.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p = inputParser;
    addRequired(p, 'tbl', @istable); % tbl with columns: conditions, spikespikes_imec0, and STIM_ON
    addParameter(p, 'PROBE_INDEX', 0, @isnumeric); % probe number
    addParameter(p, 'FIRST_BIN', 0, @isnumeric); % first time bin, aligned to STIM_ON
    addParameter(p, 'BIN_WIDTH', 50, @isnumeric); % width of each bin (ms)
    addParameter(p, 'BIN_STEP', 10, @isnumeric); % how much each bin steps/overlaps (ms)
    addParameter(p, 'N_BINS', 24, @isnumeric); % how many bins to plot
    addParameter(p, 'UNITS', [], @isnumeric);
    addParameter(p, 'IS_FHC', false, @islogical);
    addParameter(p, 'NET_THRESH', 0, @isnumeric);

    % Parse inputs
    parse(p, tbl, varargin{:});
    tbl = p.Results.tbl;
    PROBE_INDEX = p.Results.PROBE_INDEX;
    FIRST_BIN = p.Results.FIRST_BIN;
    BIN_WIDTH = p.Results.BIN_WIDTH;
    BIN_STEP = p.Results.BIN_STEP;
    N_BINS = p.Results.N_BINS;
    UNITS = p.Results.UNITS;
    IS_FHC = p.Results.IS_FHC;
    NET_THRESH = p.Results.NET_THRESH;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stimuli = cell2mat(vertcat(tbl.conditions{:}));
    xvals = sort(unique(stimuli(:,1))).';
    yvals = sort(unique(stimuli(:,2))).';

    bin_edges = arrayfun(@(x) [(FIRST_BIN + ((x*BIN_STEP)-BIN_STEP)),(FIRST_BIN + ((x*BIN_STEP)-BIN_STEP) + BIN_WIDTH)], 1:N_BINS, 'UniformOutput', false);

    prb_name = sprintf('spiketimes_%d', PROBE_INDEX);

    if isempty(UNITS)
        UNITS = 1:length(tbl.(prb_name){1});
    end

    if IS_FHC
        UNITS = 1;
    end

    all_FRs = cell(length(UNITS),1);
    for unit = 1:length(UNITS)
        FRs = cell(length(yvals),length(xvals),length(bin_edges));
        for bin = 1:length(bin_edges)
            for trial = 1:height(tbl)
                stim_ons = tbl.STIM_ON{trial};
                if IS_FHC
                    spks = tbl.(prb_name){trial};
                    nets = tbl.(sprintf('netlabels_%d',PROBE_INDEX)){trial};
                    spks = spks(nets >= NET_THRESH);
                else
                    spks = tbl.(prb_name){trial}{UNITS(unit)};
                end
                for stim = 1:length(stim_ons)
                    aligned_spks = spks-stim_ons(stim);
                    %spk_hz = sum(aligned_spks>bin_edges{bin}(1) & aligned_spks<=bin_edges{bin}(2) & nets>GAMMA)*(1000/range(bin_edges{bin}));
                    spk_hz = sum(aligned_spks>bin_edges{bin}(1) & aligned_spks<=bin_edges{bin}(2))*(1000/range(bin_edges{bin}));
                    
                    if isempty(FRs{tbl.conditions{trial}{stim}(2)==yvals,tbl.conditions{trial}{stim}(1)==xvals,bin})
                        FRs{tbl.conditions{trial}{stim}(2)==yvals,tbl.conditions{trial}{stim}(1)==xvals,bin} = spk_hz;
                    else
                        FRs{tbl.conditions{trial}{stim}(2)==yvals,tbl.conditions{trial}{stim}(1)==xvals,bin} = [FRs{tbl.conditions{trial}{stim}(2)==yvals,tbl.conditions{trial}{stim}(1)==xvals,bin}, spk_hz];
                    end
                end
            end
        end
        all_FRs{unit} = FRs;
        % fprintf(sprintf('\n----Unit %.2d complete----',UNITS(unit)))
    end
    % fprintf('\n----------------------\n')
end
