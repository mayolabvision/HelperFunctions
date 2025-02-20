function qa_pursBehavKKN_2025(data,varargin)
    % qa_pursBehavKKN_2025: Function to analyze and visualize pursuit behavior data.
    %
    % This function quickly processes eye data related to pursuit behavior
    % f1 = "pure" eye traces, split into each jump ands speed 
    %
    % Inputs:
    %   - data: (char) The name of the filepath to the raw Ripple data, with no file extention 
    %           e.g., data = '/Users/kendranoneman/OneDrive/DATA/raw/emily_walter_0347a_purs1'
    %
    % Optional Parameters (can be passed as name-value pairs):
    %   - 'PURS_PREINT': (numeric, default = 25) The time (in ms) before the pursuit target onset for plotting.
    %   - 'PURS_POSTINT': (numeric, default = 210) The time (in ms) after the pursuit target onset for plotting.
    %   - 'MS_THRESH': (numeric, default = -25, range = [-100 50]) The time (in ms) aligned to target motion onset when microsaccade must be complete
    %                  e.g., if MS_THRESH = -25 then micro-saccades cannot occur after 25 ms pre-target motion onset
    %   - 'FIG_PATH': (char or string, default = '') Path to save the generated figure. If empty, the figure is not saved.
    %
    % Outputs:
    %   - None (the function generates plots and saves them if a path is provided).
    %
    % Example usage:
    %   qa_pursBehavKKN_2025('/Users/kendranoneman/OneDrive/DATA/raw/emily_walter_0347a_purs1',...
    %                        'PURS_PREINT', 25, 'PURS_POSTINT', 210, 'MS_THRESH', -25,...
    %                        'FIG_PATH', '/Users/kendranoneman/OneDrive/DATA/processed/emily_walter_0347a/figs/purs');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create an input parser
    p = inputParser;
    addRequired(p, 'data', @(x) (ischar(x) || istable(x)));
    addParameter(p, 'PURS_PREINT', 25, @isnumeric); 
    addParameter(p, 'PURS_POSTINT', 210, @isnumeric);
    addParameter(p, 'MS_THRESH', -25, @isnumeric);
    addParameter(p, 'FIG_PATH', '', @(x) ischar(x) || isstring(x));  % Default is empty string

    % Parse the inputs
    parse(p, data, varargin{:});

    % Assign parsed values to variables
    data = p.Results.data;
    PURS_PREINT = p.Results.PURS_PREINT;
    PURS_POSTINT = p.Results.PURS_POSTINT;
    MS_THRESH = p.Results.MS_THRESH; % ms aligned to target motion onset
    FIG_PATH = p.Results.FIG_PATH;

    [nev, out_ns5, ~] = extract_nevout(data, 'SPIKE_SORT', false, 'READ_LFP', false);
    tbl = format_dataTable(nev, out_ns5, 'purs');
    [~,session,~] = fileparts(data);

    %%%%%%%%%%%%%%%%
   
    % Pure pursuit only
    f1 = figure; %('Visible','off');
    f1.Position = [100 100 1500 450*numel(unique(tbl.pursuitSpeed))];
    
    [tl,num_pure] = eyeTraces_pursSplitByConditions(tbl,unique(tbl.CROSSING_TIME(~isnan(tbl.CROSSING_TIME))),'PREINT',PURS_PREINT,'POSTINT',PURS_POSTINT,'MS_THRESH',MS_THRESH,'PURE_ONLY',true);
    
    title(tl,sprintf('%s_purs',session),'fontsize',20,'interpreter','none')
    subtitle(tl, sprintf('(# of total pure pursuit trials (w/ no MS) = %d, step-ramp duration = %d ms)',num_pure,unique(tbl.CROSSING_TIME(~isnan(tbl.CROSSING_TIME)))))
    
    if ~isempty(FIG_PATH)
        print(f1, fullfile(FIG_PATH, 'eyeTraces_pureTrials.png'), '-dpng', '-r300');
    end

end