function obj=InstallDialogdetect(obj)
    parentpanel=obj.dialogdetect;

    set(parentpanel,'Units','Pixels');
    wsize=get(parentpanel,'Position');
    set(parentpanel,'Units','normalized');
    dcolor=get(parentpanel,'Backgroundcolor');
    topDial=wsize(4)-0.3*wsize(4);
    dwidth=0.9*wsize(3);
    rightDial=0.05*wsize(3);

    selmax=4;
    selmin=1;
    strch={};
    for k=1:196
        if isempty(strch)
            strch={num2str(k)};
        else
            strch=[strch {num2str(k)}];
        end
    end
    ChSel=round(str2double(strch{1}));

    htxtlistchannels = uicontrol('Parent',parentpanel,'Style','text',...
                                 'String','channels','Units','Pixels','Position',[rightDial,topDial,dwidth,20],'BackgroundColor',dcolor);

    topDial=topDial-200;                     
    hlistchannels=uicontrol('Parent',parentpanel,'Style','listbox','max',selmax,'min',selmin,...
                            'String',strch,'Units','Pixels','Position',[rightDial,topDial,dwidth,190],'BackgroundColor',dcolor,...
                            'Callback',{@(source,event)dispchdetect_Callback(source,event,obj)});

    set([parentpanel,...
         htxtlistchannels,hlistchannels],'Units','normalized');

    function dispchdetect_Callback(source,event,obj)
        str=get(source,'string');
        val=get(source,'Value');
        ChSel=round(str2double(str(val)))';
        obj.Updatepagedetect(ChSel);
    end
end          