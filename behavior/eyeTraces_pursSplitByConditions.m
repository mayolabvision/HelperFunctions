function [tl, num_pure] = eyeTraces_pursSplitByConditions(T, rt, varargin)
    % eyeTraces_pursSplitByConditions plots radial eye velocity traces 
    % separated by pursuit speed and jump conditions, highlighting trials 
    % with and without catch-up saccades.
    %
    %%% Required Inputs: %%%
    %   T   - Table containing trial data, including eye velocity, target 
    %         motion onset, pursuit speed, jump conditions, and saccade flags.
    %   rt  - Numeric value indicating reactionTime param (ms).
    %
    %%% Optional Parameters: %%%
    %   PREINT      - Time (ms) before target motion onset to include in plot 
    %                 (default = 25).
    %   POSTINT     - Time (ms) after target motion onset to include in plot 
    %                 (default = 210).
    %   PURE_ONLY   - Logical flag indicating whether to plot only trials 
    %                 without catch-up saccades (default = true).
    %
    %%% Outputs: %%%
    %   tl       - Tiled layout handle for further customization.
    %   num_pure - Number of trials classified as "pure"
    %
    %%% Example Usage: %%%
    %   [tl, num_pure] = eyeTraces_pursSplitByConditions(T, rt, ...
    %                     'PREINT', 50, 'POSTINT', 250, 'PURE_ONLY', false);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Default values for optional parameters
    defaultPREINT = 25; 
    defaultPOSTINT = 210;
    defaultPureOnly = true;

    % Create an input parser
    p = inputParser;
    addRequired(p, 'T', @istable);
    addRequired(p, 'rt', @isnumeric); 
    addParameter(p, 'PREINT', defaultPREINT, @isnumeric); 
    addParameter(p, 'POSTINT', defaultPOSTINT, @isnumeric); 
    addParameter(p, 'PURE_ONLY', defaultPureOnly, @islogical); % PURE_ONLY must be logical

    % Parse the inputs
    parse(p, T, rt, varargin{:});

    % Assign parsed values to variables
    T = p.Results.T;
    rt = p.Results.rt;
    PREINT = p.Results.PREINT;
    POSTINT = p.Results.POSTINT;
    PURE_ONLY = p.Results.PURE_ONLY;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    angs = sort(unique(T.angle)); 
    angs = [angs; angs(1)]; % Append the first angle at the end for continuity
    
    jumps = sort(unique(T.jump));
    speeds = sort(unique(T.pursuitSpeed));
    
    % Create a tiled layout for plotting
    % 6 rows, 3 columns with compact spacing and loose padding
    tl = tiledlayout(6,3);
    tl.TileSpacing = 'compact';
    tl.Padding = 'loose';
    
    cnt = 1; % Counter for subplot tiles
    num_pure = 0; % Counter for "pure" pursuit trials
    
    for s = 1:length(speeds)
        angles_eachJump = cell(1, length(jumps));
        
        for j = 1:length(jumps)
            ax(cnt) = nexttile([2,1]); % Create subplot tile
            cnt = cnt + 1;
            
            % Extract relevant trials
            validTrials = (T.result == "CORRECT" & T.jump == jumps(j) & T.pursuitSpeed == speeds(s) & T.msFlag == 0 & ~isnan(T.csFlag));
            eyeTraces = T.eyeVel(validTrials);
            csFlags = T.csFlag(validTrials);
            angles = T.angle(validTrials);
            
            % Use the appropriate onset timing variable if it exists
            if ismember('PURSUIT_TARG', T.Properties.VariableNames)
                targetOnsets = T.PURSUIT_TARG(validTrials);
            else
                targetOnsets = T.TARG_ON(validTrials);
            end
            
            % Store angles for trials without catch-up saccades
            angles_eachJump{j} = angles(~logical(csFlags));
            
            % Compute percentage of trials removed due to saccades
            percentRemoved1 = 100 * (sum(csFlags ~= 0) / sum(validTrials));
            percentRemoved2 = 100 * (sum(csFlags == 1) / sum(validTrials));
            
            % Align eye velocity traces to target motion onset
            aligned = cellfun(@(q, w) q(:, w - PREINT:w + POSTINT), eyeTraces, num2cell(targetOnsets), 'uni', 0);
            x = (1:length(aligned{1})) - PREINT; % Time axis
            ylim([0, max(speeds) * 1.8]);
            
            % Shaded region indicating reaction time window
            yLimits = ylim;
            fill([0 rt rt 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
            
            % Plot reference speed line
            yline(speeds(s), 'r--');
            hold on;
            
            % Plot each trial's eye velocity
            for t = 1:length(eyeTraces)
                this_trl = aligned{t};
                [~, rho] = cart2pol(this_trl(1, :), this_trl(2, :));
                if csFlags(t) == 0
                    plot(x, rho, 'k-');
                    num_pure = num_pure + 1;
                else
                    if PURE_ONLY
                        continue;
                    else
                        plot(x, rho, 'r-');
                    end
                end
            end
            
            xlim([-PREINT POSTINT]);
            ylim([0, max(speeds) * 1.8]);
            
            % Add titles and labels
            if s == 1
                title(sprintf('Jump = %d', jumps(j)));
            end
            if j == 1
                ylabel({'Radial eye velocity', '(deg/s)'});
            elseif j == 3
                yl = ylabel(sprintf('Speed = %d deg/s', speeds(s)));
                yl.Rotation = -90; 
                yl.VerticalAlignment = 'top';
                yl.HorizontalAlignment = 'center';
                yl.Position = [max(xlim(ax(cnt - 1))) + 15, mean(ylim(ax(cnt - 1))), 0]; 
                yl.FontWeight = 'bold';
            end
            
            % Display percentage of trials removed
            text(ax(cnt - 1), -PREINT + 10, max(speeds) * 1.8 - 2, ...
                sprintf('%0.1f%% trials removed \n (%0.1f%% CS)', percentRemoved1, percentRemoved2), ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            
            xlabel('Time aligned to target motion onset (ms)');
            prettyFig;
        end
        
        % Plot histogram of target angles for each jump condition
        for jj = 1:length(jumps)
            ax(cnt) = nexttile([1,1]);
            cnt = cnt + 1;
            
            % Define possible directions (0:45:315 degrees)
            possibleDirections = 0:45:315;
            
            % Count occurrences of each direction
            counts = histcounts(angles_eachJump{jj}, [possibleDirections - 22.5, 337.5 + 22.5]);
            
            % Create the bar plot
            bar(possibleDirections, counts, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k');
            xticks(possibleDirections);
            xlabel('Target angle (deg)');
            title(sprintf('N = %d', length(angles_eachJump{jj})), 'FontWeight', 'normal');
            
            if jj == 1
                ylabel('Count');
            end
            
            % Reference line indicating total trials per condition
            yline(height(T(T.result=='CORRECT' & T.angle==angs(1) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s),:)),'-','Color',[0.15 0.15 0.15])
            ylim([0 height(T(T.result=='CORRECT' & T.angle==angs(1) & T.jump==jumps(j) & T.pursuitSpeed==speeds(s),:))])
            prettyFig;
        end
    end
end