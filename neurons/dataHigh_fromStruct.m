function D = dataHigh_fromStruct(exp,stimOnsets,preint,postint,epochStarts,colors_eachCondition)

    unique_conditions = unique({exp.dataMaestroPlx.condition_name}.');

    cellData = struct2cell(exp.dataMaestroPlx);
    all_units = cellData(10,:);
    
    D = struct([]);
    for t=1:length(all_units)
    
        thistrial = struct2cell(all_units{t});

        spks_binned = cell2mat(cellfun(@(q) bin_spktimes(q,0,2851,1), thistrial, 'uni', 0));
        spks_binned = spks_binned(:,(stimOnsets{t}-preint):(stimOnsets{t}+postint));

        D(t).data = spks_binned; 
        D(t).condition = exp.dataMaestroPlx(t).condition_name; 
        D(t).epochStarts = epochStarts;
        D(t).epochColors = colors_eachCondition{cellfun(@(q) isequal(exp.dataMaestroPlx(t).condition_name,q),unique_conditions,'uni',1)};
            
    
    end
    
     D(arrayfun(@(D) isempty(D.data),D)) = [];

end
