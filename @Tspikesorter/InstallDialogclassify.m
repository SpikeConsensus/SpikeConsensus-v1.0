function obj=InstallDialogclassify(obj)
    parentpanel=obj.dialogclassify;

    set(parentpanel,'Units','Pixels');
    wsize=get(parentpanel,'Position');
    set(parentpanel,'Units','normalized');
    dcolor=get(parentpanel,'Backgroundcolor');
    topDial=wsize(4)-0.1*wsize(4);
    dwidth=0.9*wsize(3);
    rightDial=0.05*wsize(3);

    selmax=4;
    selmin=1;
    strch={};
    for k=1:64
        if isempty(strch)
            strch={num2str(k)};
        else
            strch=[strch {num2str(k)}];
        end
    end
    ChSel=round(str2double(strch{1}));

    strunit={};
    for k=1:obj.Params.MCluster.nbmaxunits*3
        if isempty(strunit)
            strunit={num2str(k)};
        else
            strunit=[strunit {num2str(k)}];
        end
    end
    UnitSel=round(str2double(strunit{1}));


    Fwhitened=true;

    hFwhitened = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','Whitened','Min',0,'Max',1,'Value',Fwhitened,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fwhitened_Callback(source,event,obj)});

    topDial=topDial-100;                     
    htxtlistunits = uicontrol('Parent',parentpanel,'Style','text',...
                                 'String','Clusters','Units','Pixels','Position',[rightDial,topDial,dwidth,20],'BackgroundColor',dcolor);

    topDial=topDial-200;                         
    hlistunits=uicontrol('Parent',parentpanel,'Style','listbox','max',selmax,'min',selmin,...
                         'String',strunit,'Units','Pixels','Position',[rightDial,topDial,dwidth,190],'BackgroundColor',dcolor,...
                         'Callback',{@(source,event)dispunitclassify_Callback(source,event,obj)});                    

    set([parentpanel,hFwhitened,htxtlistunits,hlistunits],'Units','normalized');                         

    function Fwhitened_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fwhitened=true;
            obj.Updatepageclassify(UnitSel,Fwhitened);
        else
            Fwhitened=false;
            obj.Updatepageclassify(UnitSel,Fwhitened);
        end
    end

    function dispunitclassify_Callback(source,event,obj)
        str=get(source,'string');
        val=get(source,'Value');
        UnitSel=round(str2double(str(val)))';
        obj.Updatepageclassify(UnitSel,Fwhitened);
        obj.clustprofile(obj.getClusterIDMultiCh(UnitSel));
        obj.editUnitID(obj.getClusterIDMultiCh(UnitSel));                               
    end                        
end