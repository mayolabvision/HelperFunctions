function [nev,out_ns5,out_ns2] = extract_nevout(NEVPATH,varargin)
% Pulls nev and out files from raw recorded datafiles 
    %
    %%% Required inputs: %%%
    %   NEVPATH    - Path to the nev files, with or without the .nev ext
    %               (e.g. PATH = '/Users/kendranoneman/data/kendra_scrappy_0113a_rfmp2')
    %
    %%% Optional parameters: %%%
    %   SPIKE_SORT  -  True/False for spike sorting
    %                  (Default: False)
    %   GAMMA       -  minimum P(spike) value for a waveform to be classified as a
    %                  spike (between 0 and 1). If you want a more lenient sort 
    %                  (i.e. allows for more noise), then choose a smaller gamma.
    %   netFolder   -  string of the folder where the network is stored; can be
    %                  an absolute or relative path. DEFAULT: '../networks'
    %                  (e.g., '/Users/kendranoneman/Packages/nasnet/networks')
    %   READ_LFP    -  True/False for reading lfp data 
    %                  (Default: False)
    %   alignPulseEnabled - True/False 
    %
    
    %%% Outputs: %%%
    %   nev  -  3-column array containing all spikes and trial codes
    %           (nev(:,1) = channel number, nev(:,2) = sort code, nev(:,3) = time w/ 30 kHz Fs)
    %            if nev(t,1) == 0 then nev(t,2) contains trial code 
    %   out  -  structure with times, trial codes, eye data, etc...
    %           (length(dat) = number of trials)
    %
    %%% Example usage: %%%
    %   [nev,dat] = extract_nevdat('/Users/kendranoneman/data/kendra_scrappy_0113a_rfmp2'',...
    %                              'EYE_CHAN', [1,2,4], 'PUPIL_CHAN', 4, 'DIODE_CHAN', 3)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Default values for optional parameters
    defaultGamma = 0.2;
    defaultNetFolder = '../networks';

    % Create an input parser
    p = inputParser;
    addRequired(p, 'NEVPATH', @ischar); % PATH must be a string
    addParameter(p, 'SPIKE_SORT', false, @islogical); % CORRECT_ONLY must be logical
    addParameter(p, 'GAMMA', defaultGamma, @(x) isnumeric(x) && x>0 && x<1); % GAMMA must be numeric and between 0 and 1
    addParameter(p, 'netFolder', defaultNetFolder, @(x) ischar(x)); % GAMMA must be numeric and between 0 and 1
    addParameter(p, 'READ_LFP', false, @islogical); % CORRECT_ONLY must be logical
    addParameter(p, 'alignPulseEnabled', false)

    % Parse the inputs
    parse(p, NEVPATH, varargin{:});

    % Assign parsed values to variables
    NEVPATH = p.Results.NEVPATH;
    SPIKE_SORT = p.Results.SPIKE_SORT;
    GAMMA = p.Results.GAMMA;
    netFolder = p.Results.netFolder;
    READ_LFP = p.Results.READ_LFP;
    ALIGN_PULSE = p.Results.alignPulseEnabled;

    [file_path,file_name,~] = fileparts(NEVPATH);
    NEVPATH = fullfile(file_path,file_name);

    out_ns5 = read_nsx([NEVPATH '.ns5']);
    if READ_LFP
        out_ns2 = read_nsx([NEVPATH '.ns2']);
    else
        out_ns2 = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPIKE_SORT
        [slabel,nev1,net_labels,waveforms] = runNASNet([NEVPATH '.nev'],GAMMA,'netFolder',netFolder,'netname','UberNet_N50_L1');
        nev = struct();
        nev.nev = nev1;
        nev.slabel = slabel;
        nev.net_labels = net_labels;
        nev.waveforms = waveforms;
    else
        %nev = read_nev([NEVPATH '.nev']);
        nev = readNEV([NEVPATH '.nev']);
        if ~ALIGN_PULSE
            nev = nev(nev(:,2)~=0,:);
        end
    end
end