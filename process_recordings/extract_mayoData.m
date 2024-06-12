function tbl = extract_mayoData(datafolder,datafile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for extracting and preprocessing eye traces from .ns5 files
%%%%%%%%%% INPUTS %%%%%%%%%%%
% datafolder = string/char of full path where data is stored
% datafile = datafile name

%%%%%%%%%% OUTPUTS %%%%%%%%%%%
% tbl --> table with eye traces (Fs of 1000 Hz) and standard info for each trial
%    - trialName = session name + (number starting from 1 to total trials)
%    - trialOutcome = result based on eventCodes
%    - conditions (e.g. distance, angle, fixDuration, delay) = conditions separated into columns
%    - CODE TIMES (e.g. START_TRIAL, FIX_ON) = "times" in ms (but really indices), where 
%                     time = 1 is "START_TRIAL" for each trial/row 
%    - eyePos = [HE VE] eye position, after "smoothing" filtered data
%                     (filtering data with Savitzky-Golay (SG) filter)
%    - eyeVel = [HE VE] eye velocity
%    - eyeAcc = [HE VE] eye acceleration 
%    - spks   =  spike times for one channel, aligned to trial onset 

%%%%%%%%%% EXAMPLE %%%%%%%%%%%
% e.g. datafolder = '/Users/kendranoneman/Projects/mayo/HelperFunctions/process_recordings/example_data';
%      datafile = 'sb01pursA65650026';
% tbl = extract_mayoData(datafolder,datafile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% CHANGE FOR YOUR PATH %%%%%%%%%%%%%%%%%%%
% 1a. 'nevutils' package (Mayo Lab GitHub, forked from Smith Lab) should be in your MATLAB path
% https://github.com/mayolabvision/nevutils
% 1b. 'HelperFunctions' repository (which this function is contained in) should be in your MATLAB path too
% https://github.com/mayolabvision/HelperFunctions
addpath(genpath('/Users/kendranoneman/Packages')) % add nevUtils and HelperFunctions to path
addpath(genpath('/Users/kendranoneman/Projects/mayo/helperfunctions')) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Extract traces using 'nev2dat' function 
[dat,~] = nev2dat(sprintf('%s/%s',datafolder,datafile),'readNS5',true,'convertEyes',true,'include_0_255',true);
trialNames = cellfun(@(q) [datafile,char('.'),char(string(q))], num2cell(1:length(dat))','uni',0);

% 3. Pull out conditions from trial names, separated by ';' delimeter
condNames = {dat.text}';
pattern = '([^0-9;]+)(?==)';
matches = cellfun(@(q) regexp(q, pattern, 'match'), condNames, 'uni', 0);
cols = matches{1};

conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), {dat.text}.', 'uni', 0);
conditions = vertcat(conditions{:}); 

% 4. Pull eye traces, spike times, and trial codes for each trial (cell for each trial)
eye = {dat.eyedata}.'; % 3 x N (HE,VE,DI) x (N time points)
spks = cellfun(@(q,r,t) (cumsum([r; q]) - t)*1000, {dat.spiketimesdiff}.', {dat.firstspike}.', cellfun(@(q) q(1), {dat.time}.', 'uni', 0), 'uni', 0); 
trialcodes = {dat.trialcodes}.'; % C x 3 (codes) x (chan,code,time)
results = {dat.result}.';
resultNames = convertBetween_eventCodes_eventNames(results);

% 5. Make array of times using start/end time of each trial, helpful for
% aligning with trial codes and indexing eye data
times = cellfun(@(q) round(q(1)*1000:q(2)*1000), {dat.time}.', 'uni', 0);

% 6. For trial code you want to align data to (could be trial start, target
% onset, etc...), pull out eye traces around that point (preint, postint)
trialStarts = cellfun(@(q,r) find(q == round(r(r(:,2)==1,3)*1000)), times, trialcodes, 'uni', 0);

eventCodes = {dat.trialcodes}.'; eventCodes = vertcat(eventCodes{:});
eventCodes = num2cell(sort(unique(eventCodes(:,2))))';
eventCodes = eventCodes(cell2mat(eventCodes)<2000);

eventNames = convertBetween_eventCodes_eventNames(eventCodes);
trialMarkers = cellfun(@(t) cellfun(@(q,r) find(q == round(r(r(:,2)==t,3)*1000)), times, trialcodes, 'uni', 0), eventCodes, 'uni', 0)';
trialMarkers = horzcat(trialMarkers{:});
trialMarkers(cellfun('isempty',trialMarkers)) = {NaN};

% 7. Smooth the eye traces for approximating velocity and acceleration
eyePos = cellfun(@(x,y) filterEyeTraces_sglolay(x.trial(1:2,y:end),2,9,0), eye, trialStarts, 'uni', 0);
%eyePos = cellfun(@(x,y) smoothdata(x.trial(1:2,y:end),2,'gaussian',kernel), eye, trialStarts, 'uni', 0);
X = cellfun(@(q) 1:size(q,2), eyePos, 'uni', 0);

eyeVel = calcDerivative_eyeAcceleration(eyePos);
eyeAcc = calcDerivative_eyeAcceleration(eyeVel);
% 8. Save conditions and eye traces for each trial to a table

% rewrite column names based on text in dat.text
columnNames = ["trialName","trialOutcome",cols,string(eventNames)',"eyePos","eyeVel","eyeAcc","spks"];

%columnNames = ["trialName","trialOutcome","fixDuration","pursuitSpeed","angle","jumpSize",string(eventNames)',"eyePos","eyeVel","eyeAcc"];
tbl = cell2table([trialNames resultNames conditions trialMarkers cellfun(@(x){x},eyePos) cellfun(@(x){x},eyeVel) cellfun(@(x){x},eyeAcc) cellfun(@(x){x},spks)],'VariableNames',columnNames);
tbl.trialName = categorical(string(tbl.trialName)); tbl.trialOutcome = categorical(string(tbl.trialOutcome));
if ismember('emptyCnd', tbl.Properties.VariableNames)
    % Column exists, so remove it
    tbl = removevars(tbl, 'emptyCnd');
end

% 9. Save the table if you wish
save(sprintf('%s/%s_processed.mat',datafolder,datafile),'tbl','-v7');

end
