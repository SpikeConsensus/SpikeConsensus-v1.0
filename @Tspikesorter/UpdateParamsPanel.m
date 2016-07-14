function obj=UpdateParamsPanel(obj)
    if ishandle(obj.paramspanel)
        delete(obj.paramspanel);
    end
    set(obj.SpkSortoolbox,'Units','Pixels');
    wsize=get(obj.SpkSortoolbox,'Position');
    topDial=wsize(4)-0.1*wsize(4);            
    pheight=0.1*wsize(4);
    pwidth=0.995*wsize(3);
    obj.paramspanel=uipanel('Parent',obj.SpkSortoolbox,'Units','pixels','Position',[0.005*wsize(3),topDial,pwidth,pheight]);

    topDial=0.7*pheight;
    rightDial=0.005*pwidth;
    htxtRECORDING = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String','RECORDING',...
                             'Position',[rightDial,topDial,100,20],'HorizontalAlignment','left');
    topDial=0.4*pheight;
    rightDial=0.01*pwidth;
    htxtfilename = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String',['file: ' obj.filename],...
                             'Position',[rightDial,topDial,150,30],'HorizontalAlignment','left');
    topDial=0.2*pheight;            
    htxtchannelnb = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String',[num2str(obj.Params.fileinfo.Channelcount) 'channels'],...
                             'Position',[rightDial,topDial,100,30],'HorizontalAlignment','left');
    topDial=0.01*pheight;            
    htxtelectrodenb = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String',[num2str(size(obj.Params.fileinfo.electrode,2)) 'electrodes'],...
                             'Position',[rightDial,topDial,100,30],'HorizontalAlignment','left');

    topDial=0.7*pheight;
    rightDial=rightDial+pwidth/4; 
    htxtFILTERING = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String','FILTERING',...
                             'Position',[rightDial,topDial,100,20],'HorizontalAlignment','left');
    topDial=0.4*pheight;
    htxtcutoff1 = uicontrol('Parent',obj.paramspanel,'Style','text',...
                           'String',['low cutoff: ' num2str(obj.Params.detect.FcLow) 'Hz'],...
                           'Position',[rightDial,topDial,200,30],'HorizontalAlignment','left');
    topDial=0.2*pheight;                   
    htxtcutoff2 = uicontrol('Parent',obj.paramspanel,'Style','text',...
                           'String',['high cutoff: ' num2str(obj.Params.detect.FcHigh) 'Hz'],...
                           'Position',[rightDial,topDial,200,30],'HorizontalAlignment','left');
    topDial=0.01*pheight;                   
    htxtcutoff3 = uicontrol('Parent',obj.paramspanel,'Style','text',...
                           'String',['Sample freq: ' num2str(obj.Params.fileinfo.samplingrate/obj.dwsamplefactor) 'Hz'],...
                           'Position',[rightDial,topDial,200,30],'HorizontalAlignment','left');

    topDial=0.7*pheight;
    rightDial=rightDial+pwidth/4;                    
    htxtDETECTION = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String','DETECTION',...
                             'Position',[rightDial,topDial,100,20],'HorizontalAlignment','left');                   
    topDial=0.2*pheight;                    
    htxtdetect = uicontrol('Parent',obj.paramspanel,'Style','text',...
                           'String',['z-threshold: ' num2str(obj.Params.detect.zthresh) 'std'],...
                           'Position',[rightDial,topDial,200,30],'HorizontalAlignment','left');
    topDial=0.7*pheight;
    rightDial=rightDial+pwidth/4;                    
    htxtCLUSTERING = uicontrol('Parent',obj.paramspanel,'Style','text',...
                             'String','CLUSTERING',...
                             'Position',[rightDial,topDial,100,20],'HorizontalAlignment','left');                                  
    set([obj.SpkSortoolbox,obj.paramspanel,...
         htxtRECORDING,htxtfilename,htxtchannelnb,htxtelectrodenb,...
         htxtFILTERING,htxtcutoff1,htxtcutoff2,htxtcutoff3,...
         htxtDETECTION,htxtdetect,...
         htxtCLUSTERING],'Units','normalized'); 
end