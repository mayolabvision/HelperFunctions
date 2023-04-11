function sliderChangingFcn(~,event,data1,data2,Obj1a,Obj1b,Obj2a,Obj2b,Ell1a,Ell1b,Ell2a,Ell2b,final_table,colors)
    % Update heatmap and title with new selection of data
    value = round(event.Value);
    ObjA = polarplot(data1(final_table.Type=="pure",value,1),data1(final_table.Type=="pure",value,2), ...
                     'o','color',colors{3},'markersize',4,'linewidth',2,'visible','off','handlevisibility','off');
    Obj1a.ThetaData  =  ObjA.ThetaData; Obj1a.RData      =  ObjA.RData;
    [x1,y1] = pol2cart(data1(final_table.Type=="pure",value,1),data1(final_table.Type=="pure",value,2)); 
    [r_ellipse1, X0a, Y0a] = error_ellipseJPM([x1,y1]);
    [t1,r1] = cart2pol(r_ellipse1(:,1)+X0a,r_ellipse1(:,2)+Y0a);
    EllA = polarplot(t1,r1,'color',colors{3},'linewidth',2,'visible','off','handlevisibility','off');
    Ell1a.ThetaData  =  EllA.ThetaData; Ell1a.RData      =  EllA.RData;
    
    ObjB = polarplot(data1(final_table.Type=="forward",value,1),data1(final_table.Type=="forward",value,2), ...
                    'o','color',colors{2},'markersize',4,'linewidth',2,'visible','off','handlevisibility','off');
    Obj1b.ThetaData  =  ObjB.ThetaData; Obj1b.RData      =  ObjB.RData;
    [x1,y1] = pol2cart(data1(final_table.Type=="forward",value,1),data1(final_table.Type=="forward",value,2)); 
    [r_ellipse1, X0a, Y0a] = error_ellipseJPM([x1,y1]);
    [t1,r1] = cart2pol(r_ellipse1(:,1)+X0a,r_ellipse1(:,2)+Y0a);
    EllB = polarplot(t1,r1,'color',colors{2},'linewidth',2,'visible','off','handlevisibility','off');
    Ell1b.ThetaData  =  EllB.ThetaData; Ell1b.RData      =  EllB.RData;
    

    % velocity
    ObjA = polarplot(data2(final_table.Type=="pure",value,1),data2(final_table.Type=="pure",value,2), ...
                     'o','color',colors{3},'markersize',4,'linewidth',2,'visible','off','handlevisibility','off');
    Obj2a.ThetaData  =  ObjA.ThetaData; Obj2a.RData      =  ObjA.RData;
    [x1,y1] = pol2cart(data2(final_table.Type=="pure",value,1),data2(final_table.Type=="pure",value,2)); 
    [r_ellipse1, X0a, Y0a] = error_ellipseJPM([x1,y1]);
    [t1,r1] = cart2pol(r_ellipse1(:,1)+X0a,r_ellipse1(:,2)+Y0a);
    EllA = polarplot(t1,r1,'color',colors{3},'linewidth',2,'visible','off','handlevisibility','off');
    Ell2a.ThetaData  =  EllA.ThetaData; Ell2a.RData      =  EllA.RData;
    

    ObjB = polarplot(data2(final_table.Type=="forward",value,1),data2(final_table.Type=="forward",value,2), ...
                    'o','color',colors{2},'markersize',4,'linewidth',2,'visible','off','handlevisibility','off');
    Obj2b.ThetaData  =  ObjB.ThetaData; Obj2b.RData      =  ObjB.RData;
    [x1,y1] = pol2cart(data2(final_table.Type=="forward",value,1),data2(final_table.Type=="forward",value,2)); 
    [r_ellipse1, X0a, Y0a] = error_ellipseJPM([x1,y1]);
    [t1,r1] = cart2pol(r_ellipse1(:,1)+X0a,r_ellipse1(:,2)+Y0a);
    EllB = polarplot(t1,r1,'color',colors{2},'linewidth',2,'visible','off','handlevisibility','off');
    Ell2b.ThetaData  =  EllB.ThetaData; Ell2b.RData      =  EllB.RData;
    
end