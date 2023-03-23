function [exp_clean,unique_conditions] = struct_pullConditions(exp,extract_conditions,extract_columns,define_columns)
% OBJECTIVE:
% pull out trials with particular conditions, based on trType field
%
% INPUTS:
% exp = where data in rows of exp.dataMaestroPlx are individual trials
% extract_conditions = the names of conditions you want to keep, in columns where you want to exclude some of the conditions
%                      e.g. {'1fXXX','2fXXX'}
% extract_columns = indices of columns that "extract_conditions" come from
%                      e.g. [3 4]
% define_columns = indices of columns you want to label trials with 
%                      e.g. [1 2 5]
%
% OUTPUTS:
% exp_clean = new exp including trials w/ specified conditions
% unique_conditions = based on the "define columns", the conditions you have
%                     e.g. {'d000';'d090';'d180';'d270'}

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

if nargin < 4
    define_columns = 1:size(trTypes,1);
end 

conditions = trTypes(:,define_columns);
combine_conditions = cell(size(conditions,1),1);
for c=1:size(conditions,2)
    if c==1
        combine_conditions = cellfun(@(x,y)[x,y], combine_conditions,conditions(:,c),'uni',0); 
    else
        combine_conditions = cellfun(@(x,y)[x,'_',y], combine_conditions,conditions(:,c),'uni',0); 
    end
end

[exp_clean.dataMaestroPlx.condition_name] = combine_conditions{:};
unique_conditions = unique(combine_conditions);

end

