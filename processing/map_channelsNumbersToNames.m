function [mappings,probe_specs] = map_channelsNumbersToNames(mapping_name,probeID,varargin)
    % dfds
    % kjlad
    function [rel_depths,abs_depths] = getChannelDepths(chan_num,distanceToTip,interElecSpacing,probeDepth)
        rel_depths = zeros(length(chan_num),1);
        num_chans = length(chan_num);
        for c = 1:num_chans
            rel_depths(c) = (0 - interElecSpacing*((num_chans-c)) - distanceToTip)./1000;
        end
        abs_depths = rel_depths + probeDepth;
    end

    % Default values for optional parameters
    [currentDir, ~, ~] = fileparts(mfilename('fullpath'));
    defaultMappingPath = fullfile(fileparts(currentDir), 'mappings');

    % Create an input parser
    p = inputParser;
    addRequired(p, 'mapping_name', @(x) (iscategorical(x) ||  isstring(x) || ischar(x)));
    addRequired(p, 'probeID', @(x) (iscell(x) ||  isstring(x) || ischar(x) || iscategorical(x)));
    addParameter(p, 'mappingPath', defaultMappingPath, (@(x) ischar(x))); 
    addParameter(p, 'probeDepths_mm', 0, (@(x) (isinteger(x) || iscell(x))));

    % Parse the inputs
    parse(p, mapping_name, probeID, varargin{:});

    % Assign parsed values to variables
    mapping_name = p.Results.mapping_name;
    probeID = p.Results.probeID;
    mappingPath = p.Results.mappingPath;
    probeDepths_mm = p.Results.probeDepths_mm;

    if ~isequal(class(probeDepths_mm),'cell')
        if probeDepths_mm==0
            probeDepths_mm = num2cell(repmat(probeDepths_mm,1,length(probeID)+1));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iscategorical(mapping_name)
        mapping_name = char(string(mapping_name));
    end

    mappings = readtable([mappingPath,'/',mapping_name,'.csv']);
    mappings = sortrows(mappings, 'chan_num');
    mappings.mapped_name = categorical(mappings.mapped_name);
    
    probes = readtable([mappingPath,'/','probe_configs.csv']);

    probe_specs = [];
    for probe = 1:length(probeID)
        this_probe = probes(ismember(probes.probeID,probeID{probe}),:);
        probe_specs = [probe_specs; this_probe];
        [rel_depths, abs_depths] = getChannelDepths(mappings.chan_num,this_probe.distanceTipToFirstElectrode,this_probe.interElectrodeSpacing,probeDepths_mm{probe}); 
    end
    probe_specs = table2struct(probe_specs);
    probe_specs.numberUses = [];

    mappings.relDepth_mm = rel_depths;
    mappings.absDepth_mm = abs_depths;

end