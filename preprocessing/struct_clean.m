function [exp_clean,unitnames,snrs] = struct_clean(exp)
% OBJECTIVE:
% clean up struct to not include empty trials or bad units, and reorder unit names to consistent across trials
%
% INPUTS:
% exp = where data in rows of exp.dataMaestroPlx are individual trials
%
% OUTPUTS:
% exp_clean = new exp including trials w/ specified conditions

exp_clean = exp;

% throw out empty trials
exp_clean.dataMaestroPlx(find(cellfun(@isempty, {exp_clean.dataMaestroPlx.units}.'))) = [];

% throw out trials missing stimulus onset time
tagS  =  {exp_clean.dataMaestroPlx.tagSection}.'; tagS = vertcat(tagS{:});
if sum(cellfun(@(q) isempty(q), {tagS.stTimeMS}.', 'uni', 1))~=length(tagS) % check that field exists
    exp_clean.dataMaestroPlx(cellfun(@(q) isempty(q), {tagS.stTimeMS}.', 'uni', 1)) = [];
end

% find names of units that don't drop over course of session
[unitnames,snrs] = findConsistentUnits_fromStruct(exp_clean);
exp_clean.info.channels = unitnames; exp_clean.info.SNRs = snrs; % replace channels/snrs with new names/order

% toss out "bad" units and sort units in numerical/alphabetical order
spk_cnts = zeros(length(exp_clean.dataMaestroPlx),1);
for t=1:length(exp_clean.dataMaestroPlx)
    exp_clean.dataMaestroPlx(t).units  =  rmfield(exp_clean.dataMaestroPlx(t).units,setdiff(fieldnames(exp_clean.dataMaestroPlx(t).units),unitnames));
    exp_clean.dataMaestroPlx(t).units  =  orderfields(exp_clean.dataMaestroPlx(t).units,unitnames);

    sptimes = struct2cell(exp_clean.dataMaestroPlx(t).units);
    spk_cnts(t) = mean(cellfun(@length, sptimes));

    % replace condition w/ "rotated" direction
    newdir = wrapTo360(str2double(exp_clean.dataMaestroPlx(t).trType(strfind(exp_clean.dataMaestroPlx(t).trType,'d')+1:strfind(exp_clean.dataMaestroPlx(t).trType,'d')+3))-(-1*double(exp_clean.info.rotfactor)));
    exp_clean.dataMaestroPlx(t).trType = strrep(exp_clean.dataMaestroPlx(t).trType,exp_clean.dataMaestroPlx(t).trType(strfind(exp_clean.dataMaestroPlx(t).trType,'d')+1:strfind(exp_clean.dataMaestroPlx(t).trType,'d')+3),sprintf('%03d',newdir));
end

% remove trials where the mean spike count exceeded 3 standard deviations from mean spike count
exp_clean.dataMaestroPlx(spk_cnts < (mean(spk_cnts)-3*std(spk_cnts)) | spk_cnts > (mean(spk_cnts)+3*std(spk_cnts))) = [];

end