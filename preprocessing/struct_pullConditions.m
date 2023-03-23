function [exp_clean,condition_names] = struct_pullConditions(exp,extract_conditions,extract_columns,define_columns)
% OBJECTIVE:
% determine if nhp made microsaccade around stimulus onset in given trial
%
% INPUTS:
% eye = 4x1 cell array w/ (HEPos, VEPos, HEVel, VEVel)
% stimOnset = time of target motion onset or other stimulus onset
% preint = time before stimOnset you want to include 
% postint = time after stimOnset you want to include
% accThresh = acceleration threshold for microsacc detection (e.g. 750 deg/s^2)
% velThresh = velocity threshold for microsacc detection (e.g. 50 deg/s)
%
% OUTPUTS:
% msFlag = 1 if microsacc detected, 0 if not detected

% if nargin < 3
%     preint = 50; postint = 50; accThresh = 750; velThresh = 50;
% elseif nargin < 5
%     accThresh = 750; velThresh = 50;
% end
exp_clean = exp;

% pull out conditions of interest
trTypes_all  =  {exp_clean.dataMaestroPlx.trType}.';
trTypes  =  cellfun(@(x) cellstr(strsplit(x, '_')), trTypes_all, 'uni', 0);
trTypes  =  vertcat(trTypes{:});
trials_include = sum(ismember(trTypes,extract_conditions),2) == length(extract_columns);

exp_clean.dataMaestroPlx(~trials_include) = [];

% rename conditions based on values you care about
trTypes_all  =  {exp_clean.dataMaestroPlx.trType}.';
trTypes  =  cellfun(@(x) cellstr(strsplit(x, '_')), trTypes_all, 'uni', 0);
trTypes  =  vertcat(trTypes{:});

conditions = trTypes(:,define_columns);

if length(define_columns)==1
    [exp_clean.dataMaestroPlx.condition_name] = conditions{:};
else
    
end

if length(conditionName_columns)==1
    conditions = trTypes(:,conditionName_columns);
elseif length(conditionName_columns)==2
    conditions = cellfun(@(x,y)[x,'_',y], trTypes(:,conditionName_columns(1)),trTypes(:,conditionName_columns(2)),'uni',1);
elseif length(conditionName_columns)==3
    conditions = cellfun(@(x,y,z)[x,'_',y,'_',z], trTypes(:,conditionName_columns(1)),trTypes(:,conditionName_columns(2)),trTypes(:,conditionName_columns(3)),'uni',1);
end
    
unique_conditions = unique(conditions);

end

