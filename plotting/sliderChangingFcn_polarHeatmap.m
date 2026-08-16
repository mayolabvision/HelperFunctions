function sliderChangingFcn_polarHeatmap(~,event,data,heatObj1,heatObj2,final_table)
    % Update heatmap and title with new selection of data
    value = round(event.Value);
    heatObjA = histogram2Polar(data(final_table.Type=="pure",value,1), data(final_table.Type=="pure",value,2), 0.1, ...
                               'MapColor',flip(gray),'ThetaZeroLocation','right','RLim',[0 2]);
    heatObj1.Data = heatObjA.Data;

    heatObjB = histogram2Polar(data(final_table.Type=="forward",value,1), data(final_table.Type=="forward",value,2), 0.1, ...
                               'MapColor',flip(summer),'ThetaZeroLocation','right','RLim',[0 2]);
    heatObj2.Data = heatObjB.Data;
end