clear
clc

tic
%% Setup folders and parameters
dataFolder  =  '/Users/kendranoneman/Projects/mayo/data/mt-motion_pulse/raw';
figFolder   =  '/Users/kendranoneman/Projects/mayo/figs/mt-motion_pulse';
tableFolder  =  '/Users/kendranoneman/Projects/mayo/data/mt-motion_pulse/processed';

addpath(genpath('/Users/kendranoneman/Projects/mayo/HelperFxns'))

filelist   =  dir(dataFolder); % files in raw data folder
name       =  {filelist.name}; 
sessions   =  name(~strncmp(name, '.', 1))';

preint   =  300; % time before onset
postint  =  600; % time after onset 
motionStart  =  1250; % constant across trials

types        =  {'pure','forward','backward'}; % types of trials 

%%
parpool('local',8);
warning('off')
[final_table] = cell(length(sessions),5);
parfor s=1:length(sessions) % loop through each session
    sess_name  =  string(sessions{s}); sess_name = char(extractBetween(sess_name,".","."));
    fprintf(sprintf('%d / 37 \n\n',s));

    % Monkey name
    if isequal(sess_name(2),'a')
        monk  =  'aristotle';
    elseif isequal(sess_name(2),'b')
        monk  =  'batman';
    end

    dat = load('-mat',sprintf('%s/%s',dataFolder,sessions{s})); % raw structure
    [exp_clean,unitnames,snrs] = struct_clean(dat.exp);

    tagS  =  {exp_clean.dataMaestroPlx.tagSection}.'; tagS = vertcat(tagS{:});
    if sum(cellfun(@(q) isempty(q), {tagS.stTimeMS}.', 'uni', 1))~=length(tagS)
        stmFlag = 1;
        stimOnsets = {tagS.stTimeMS}.';
    else
        stmFlag = 0;
        stimOnsets = num2cell(repmat(motionStart,length(exp_clean.dataMaestroPlx),1));
    end

    % clean up struct
    [msFlag,eye_adjust] = cellfun(@(q,m) detect_msTrials(struct2cell(q),m,50,100,750,50), {exp_clean.dataMaestroPlx.mstEye}.', stimOnsets, 'uni', 0);
    exp_clean.dataMaestroPlx(logical(cell2mat(msFlag))) = []; eye_adjust = eye_adjust(~logical(cell2mat(msFlag))); 
    [exp_clean.dataMaestroPlx.mstEye] = eye_adjust{:};

    extract_conditions = {'1fXXX','2fXXX'};
    extract_columns = [3 4];
    define_columns = 1;
    [exp_clean,condition_names] = struct_pullConditions(exp_clean,extract_conditions,extract_columns,define_columns);

    %%%%%%%%%%%%%%%%% TRIAL TABLE %%%%%%%%%%%%%%%%%
    if stmFlag==1
        tagS  =  {exp_clean.dataMaestroPlx.tagSection}.'; tagS = vertcat(tagS{:});
        stimOnsets = {tagS.stTimeMS}.';
    else
        stimOnsets = num2cell(repmat(motionStart,length(exp_clean.dataMaestroPlx),1));
    end
    [pursuitOnsets,rxnTimes] = cellfun(@(q,m) detect_pursuitOnset(struct2cell(q),m,50,300), {exp_clean.dataMaestroPlx.mstEye}.', stimOnsets, 'uni', 0);
    exp_clean.dataMaestroPlx(isnan(cell2mat(rxnTimes))) = []; pursuitOnsets(isnan(cell2mat(rxnTimes))) = []; stimOnsets(isnan(cell2mat(rxnTimes))) = []; rxnTimes(isnan(cell2mat(rxnTimes))) = []; 

    motionDirs = cellfun(@(q) str2double(q(2:4)), {exp_clean.dataMaestroPlx.trType}.', 'uni', 0);
    [csTypes,ipt,saccProps] = cellfun(@(q,p,d) detect_catchupSaccade(struct2cell(q),p,d,1,200,750,30), {exp_clean.dataMaestroPlx.mstEye}.', pursuitOnsets, motionDirs, 'uni', 0);
    exp_clean.dataMaestroPlx(cell2mat(csTypes)==0) = []; pursuitOnsets(cell2mat(csTypes)==0) = []; stimOnsets(cell2mat(csTypes)==0) = []; rxnTimes(cell2mat(csTypes)==0) = []; motionDirs(cell2mat(csTypes)==0) = []; ipt(cell2mat(csTypes)==0) = []; saccProps(cell2mat(csTypes)==0) = []; csTypes(cell2mat(csTypes)==0) = [];
    new_condition = cellfun(@(x,y)[x,'_',types{y}], {exp_clean.dataMaestroPlx.condition_name}.',csTypes,'uni',0);
    [exp_clean.dataMaestroPlx.condition_name] = new_condition{:};

    tt = [cellstr(repmat(categorical(string(monk)),length(rxnTimes),1)) cellstr(repmat(categorical(string(sess_name)),length(rxnTimes),1)) ...
                   {exp_clean.dataMaestroPlx.trName}.' {exp_clean.dataMaestroPlx.trType}.' {exp_clean.dataMaestroPlx.condition_name}.'...
                   motionDirs cellfun(@(q) types{q}, csTypes, 'uni', 0) stimOnsets pursuitOnsets rxnTimes ipt saccProps ...
                   cellfun(@(q) struct2cell(q), {exp_clean.dataMaestroPlx.mstEye}.','uni',0) ...
                   cellfun(@(r) r(1:2), cellfun(@(q) struct2cell(q), {exp_clean.dataMaestroPlx.target}.', 'uni', 0), 'uni', 0)];

    varNames   =  ["Monkey","Session","TrialName","TrialType","TrialCondition","Direction","Type","TargetMotionOnset","PursuitOnset","RxnTime","SaccadeTimes","SaccadeProps","EyeTraces","TargetTraces"];
    trialTbl   =  cell2table(tt,"VariableNames",varNames); 
    trialTbl.Monkey     =  categorical(string(trialTbl.Monkey)); trialTbl.Session    =  categorical(string(trialTbl.Session));
    trialTbl.TrialName   =  categorical(string(trialTbl.TrialName)); trialTbl.TrialType   =  categorical(string(trialTbl.TrialType));
    trialTbl.TrialCondition   =  categorical(string(trialTbl.TrialCondition)); trialTbl.Type   =  categorical(string(trialTbl.Type));
    %%%%%%%%%%%%%%%%%% DataHigh %%%%%%%%%%%%%%%%%
    % 0_b, 0_f, 0_p, 90_b, 90_f, 90_p, etc...
    colors_eachCondition = {[232,125,12]./255; [22,95,19]./255; [95,23,137]./255; [239,151,42]./255; [75,153,86]./255; [130,83,174]./255; [241,164,68]./255; [104,182,122]./255; [157,129,203]./255; [252,220,129]./255; [111,210,142]./255; [175,166,231]./255};
    epochStarts = 1;
    D = dataHigh_fromStruct(exp_clean,stimOnsets,preint,postint,epochStarts,colors_eachCondition);

    %%%%%%%%%%%%%%% UNIT TABLE %%%%%%%%%%%%%%%%
    ut = makeUnitsTable_fromStruct(exp_clean,trialTbl,stimOnsets,100);
    ut = [cellstr(repmat(categorical(string(monk)),length(unitnames),1)) cellstr(repmat(categorical(string(sess_name)),length(unitnames),1)) ut];
    
    varNames   =  ["Monkey","Session","UnitName","Sess_Unit","BrainArea","SNR","BestDir","NullDir","PrefDirFit","DepthMod","SelDir","DI","SI","signiffl","SpikeTimes"];
    unitTbl   =  cell2table(ut,"VariableNames",varNames); 
    unitTbl.Monkey = categorical(string(unitTbl.Monkey)); unitTbl.Session = categorical(string(unitTbl.Session));
    unitTbl.UnitName = categorical(string(unitTbl.UnitName)); unitTbl.Sess_Unit = categorical(string(unitTbl.Sess_Unit)); unitTbl.BrainArea = categorical(string(unitTbl.BrainArea));
    
    final_table(s,:) = {string(monk) string(sess_name) trialTbl, unitTbl, D};  
end
delete(gcp('nocreate'))

varNames = ["Monkey","Session","Trials","Units","DataHigh"];
final_table = cell2table(final_table,"VariableNames",varNames);
final_table.Monkey = categorical(final_table.Monkey); final_table.Session = categorical(final_table.Session);

save(sprintf('%s/ftbl-n3.mat',tableFolder),'final_table','-v7.3');

toc
load gong.mat
sound(y)
