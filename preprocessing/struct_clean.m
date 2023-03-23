function exp_clean = struct_clean(exp)

exp_clean = exp;

% throw out empty trials
exp_clean.dataMaestroPlx(find(cellfun(@isempty, {exp_clean.dataMaestroPlx.units}.'))) = [];

% throw out trials missing stimulus onset time
tagS  =  {exp_clean.dataMaestroPlx.tagSection}.'; tagS = vertcat(tagS{:});
if sum(cellfun(@(q) isempty(q), {tagS.stTimeMS}.', 'uni', 1))~=length(tagS) % check that field exists
    exp_clean.dataMaestroPlx(cellfun(@(q) isempty(q), {tagS.stTimeMS}.', 'uni', 1)) = [];
end

% find names of units that don't drop over course of session
channels       =  exp_clean.info.channels; % names of all channels
snrs           =  exp_clean.info.SNRs; % SNR for each channel 

all_units      =  cellfun(@(x) fieldnames(x), {exp_clean.dataMaestroPlx.units}.', 'uni', 0);
[B,BG]         =  groupcounts(vertcat(all_units{:}));
[~,ia]         =  setdiff(channels,cellfun(@(y) y(end-3:end), BG(B==max(B)), 'uni', 0));
channels(ia)   =  []; snrs(ia) = [];

[unitnames,I]  =  sort(channels); snrs = snrs(I);
exp_clean.info.channels = unitnames; exp_clean.info.SNRs = snrs; % replace channels/snrs with new names/order

unitnames      =  cellfun(@(z) strcat('unit',z), unitnames, 'uni', 0)';

% toss out "bad" units and sort units in numerical/alphabetical order
spk_cnts = zeros(length(exp_clean.dataMaestroPlx),1);
for t=1:length(exp_clean.dataMaestroPlx)
    exp_clean.dataMaestroPlx(t).units  =  rmfield(exp_clean.dataMaestroPlx(t).units,setdiff(fieldnames(exp_clean.dataMaestroPlx(t).units),unitnames));
    exp_clean.dataMaestroPlx(t).units  =  orderfields(exp_clean.dataMaestroPlx(t).units,unitnames);

    sptimes = struct2cell(exp_clean.dataMaestroPlx(t).units);
    spk_cnts(t) = mean(cellfun(@length, sptimes));
end

% remove trials where the mean spike count exceeded 3 standard deviations from mean spike count
exp_clean.dataMaestroPlx(spk_cnts < (mean(spk_cnts)-3*std(spk_cnts)) | spk_cnts > (mean(spk_cnts)+3*std(spk_cnts))) = [];

end