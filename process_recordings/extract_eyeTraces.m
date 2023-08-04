function tbl = extract_eyeTraces(datafolder,datafile,eventNames,kernel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for extracting and preprocessing eye traces from .ns5 files
%%%%%%%%%% INPUTS %%%%%%%%%%%
% datafolder = string/char of full path where data is stored
% datafile = datafile name
% eventNames = cell array of eventNames you want to store the times of
% kernel = how much to smooth eye position before calculating vel/acc

%%%%%%%%%% OUTPUTS %%%%%%%%%%%
% tbl --> table with eye traces (Fs of 1000 Hz) and standard info for each trial
%    - trialName = session name + (number starting from 1 to total trials)
%    - trialOutcome = result based on eventCodes
%    - conditions (e.g. fixDuration) = conditions separated into columns
%    - times (e.g. FIX_ON) = "times" in ms (but really indices), where 
%                     time = 1 is "START_TRIAL", for desired events 
%    - eyePos = [HE VE] eye position, after "smoothing" filtered data
%                     (filtering data with a Gaussian window, based on kernel input)
%    - eyeVel = [HE VE] eye velocity
%    - eyeAcc = [HE VE] eye acceleration 

%%%%%%%%%% EXAMPLE %%%%%%%%%%%
% e.g. datafolder = '/Users/kendranoneman/Projects/mayo/HelperFxns/process_recordings/example_data';
%      datafile = 'sb01pursA65650026';
%      eventNames = {'START_TRIAL','FIX_ON','FIX_OFF','REWARD','TARG_ON','FIXATE','IGNORED','BROKE_FIX','BROKE_TARG'}
% "extract_eyeTraces(dataFolder,datafile,eventNames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. 'nevutils' package (downloaded from Mayo Lab GitHub, forked from Smith
% Lab) should be already in your MATLAB path

% 2. Extract traces using 'nev2dat' function 
[dat,~] = nev2dat(sprintf('%s/%s',datafolder,datafile),'readNS5',true,'convertEyes',true);
trialNames = cellfun(@(q) [datafile,char('.'),char(string(q))], num2cell(1:length(dat))','uni',0);

% 4. Pull out conditions from trial names, separated by ';' delimeter
conditions = cellfun(@(x) cellfun(@(q,r) str2double(x(q+1:r-1)), num2cell(strfind(x,'=')), num2cell(strfind(x,';')), 'uni', 0), {dat.text}.', 'uni', 0);
conditions = vertcat(conditions{:}); conditions = conditions(:,2:5);

% 5. Pull eye traces and trial codes for each trial (cell for each trial)
eye = {dat.eyedata}.'; % 3 x N (HE,VE,DI) x (N time points)
trialcodes = {dat.trialcodes}.'; % C x 3 (codes) x (chan,code,time)
results = {dat.result}.';
results(cellfun(@isnan, results)) = {150}; % temporary fix for no correct
resultNames = convertBetween_eventCodes_eventNames(results);

% 6. Make array of times using start/end time of each trial, helpful for
% aligning with trial codes and indexing eye data
times = cellfun(@(q) round(q(1)*1000:q(2)*1000), {dat.time}.', 'uni', 0);

% 7. For trial code you want to align data to (could be trial start, target
% onset, etc...), pull out eye traces around that point (preint, postint)
trialStarts = cellfun(@(q,r) find(q == round(r(r(:,2)==1,3)*1000)), times, trialcodes, 'uni', 0);

eventCodes = convertBetween_eventCodes_eventNames(eventNames);
trialMarkers = cellfun(@(t) cellfun(@(q,r) find(q == round(r(r(:,2)==t,3)*1000)), times, trialcodes, 'uni', 0), eventCodes, 'uni', 0)';
trialMarkers = horzcat(trialMarkers{:});
trialMarkers(cellfun('isempty',trialMarkers)) = {NaN};

eyePos = cellfun(@(x,y) smoothdata(x.trial(1:2,y:end),2,'gaussian',kernel), eye, trialStarts, 'uni', 0);
X = cellfun(@(q) 1:size(q,2), eyePos, 'uni', 0);

% 9. Smooth the eye traces for approximating velocity and acceleration
eyeVel = cellfun(@(q,x) [gradient(q(1,:)') ./ gradient(x(:)./1000),  gradient(q(2,:)') ./ gradient(x(:)./1000)]', eyePos, X, 'uni', 0);
eyeAcc = cellfun(@(q,x) [gradient(q(1,:)') ./ gradient(x(:)./1000),  gradient(q(2,:)') ./ gradient(x(:)./1000)]', eyeVel, X, 'uni', 0);

% 10. Save conditions and eye traces for each trial to a table
columnNames = ["trialName","trialOutcome","fixDuration","pursuitSpeed","angle","jumpSize",string(eventNames),"eyePos","eyeVel","eyeAcc"];
tbl = cell2table([trialNames resultNames conditions trialMarkers cellfun(@(x){x},eyePos) cellfun(@(x){x},eyeVel) cellfun(@(x){x},eyeAcc)],'VariableNames',columnNames);
tbl.trialName = categorical(string(tbl.trialName)); tbl.trialOutcome = categorical(string(tbl.trialOutcome));

end