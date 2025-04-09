function qa_pursBehavKKN_2025(nevpath,varargin)
    % qa_pursBehavKKN_2025: Function to analyze and visualize pursuit behavior data.
    %
    % This function quickly processes eye data related to pursuit behavior
    % f1 = "pure" eye traces, split into each jump ands speed 
    %
    % Inputs:
    %   - nevpath: (char) The full filepath to the Ripple files, with or without .nev
    %           e.g., nevpath = '/Users/kendranoneman/OneDrive/DATA/raw/emily_walter_0347a_purs1'
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
    %   qa_pursBehavKKN_2025('/Volumes/lab_NHPdata/kendra_scrappy_0121a_purs2.nev',...
    %                        'PURS_PREINT', 25, 'PURS_POSTINT', 210, 'MS_THRESH', -25,...
    %                        'FIG_PATH','/Users/kendranoneman/OneDrive/DATA/kendra_scrappy_0121a')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an input parser
    p = inputParser;
    addRequired(p, 'nevpath', @(x) (ischar(x) || istable(x)));
    addParameter(p, 'PURS_PREINT', 25, @isnumeric); 
    addParameter(p, 'PURS_POSTINT', 210, @isnumeric);
    addParameter(p, 'MS_THRESH', -25, @isnumeric);
    addParameter(p, 'FIG_PATH', fileparts(nevpath), @(x) ischar(x) || isstring(x));  % Default is empty string

    % Parse the inputs
    parse(p, nevpath, varargin{:});

    % Assign parsed values to variables
    nevpath = p.Results.nevpath;
    PURS_PREINT = p.Results.PURS_PREINT;
    PURS_POSTINT = p.Results.PURS_POSTINT;
    MS_THRESH = p.Results.MS_THRESH; % ms aligned to target motion onset
    FIG_PATH = p.Results.FIG_PATH;

    [~,session,~] = fileparts(nevpath);

    % Processing Ripple data
    tic
    fprintf('\n---- generating nev_out for %s ----\n', session);
    [nev, out_ns5, ~] = extract_nevout(nevpath);

    fprintf('\n---- generating dat for %s ----\n', session);
    [dat,~] = format_datTrials(nev, out_ns5);

    fprintf('\n---- generating tbl for %s ----\n', session);
    tbl = convert_smithDat_mayoTbl(dat,'TASK_NAME','purs');

    %%%%%%%%%%%%%%%%
    if ~exist(FIG_PATH, 'dir'), mkdir(FIG_PATH); end
    fprintf('\n---- making purs fig for %s ----\n', session);

    % Pure pursuit only
    f1 = figure; %('Visible','off');
    f1.Position = [100 100 1500 450*numel(unique(tbl.pursuitSpeed))];
    
    [tl,num_pure] = eyeTraces_pursSplitByConditions(tbl,unique(tbl.CROSSING_TIME(~isnan(tbl.CROSSING_TIME))),'PREINT',PURS_PREINT,'POSTINT',PURS_POSTINT,'MS_THRESH',MS_THRESH,'PURE_ONLY',true);
    
    title(tl,sprintf('%s_purs',session),'fontsize',20,'interpreter','none')
    subtitle(tl, sprintf('(# of total pure pursuit trials (w/ no MS) = %d, step-ramp duration = %d ms)',num_pure,unique(tbl.CROSSING_TIME(~isnan(tbl.CROSSING_TIME)))))
    
    if ~isempty(FIG_PATH)
        print(f1, fullfile(FIG_PATH, 'eyeTraces_pureTrials.png'), '-dpng', '-r300');
    end

    tc = toc;
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(sprintf('Total elapsed time was %2.2f minutes',tc/60))
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');

end