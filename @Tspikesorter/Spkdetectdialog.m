function obj=Spkdetectdialog(obj)
    if ~ischar(obj.SpkSortoolbox)
        set(obj.SpkSortoolbox,'Units','Pixels');
        toolboxpos=get(obj.SpkSortoolbox,'Position');
        set(obj.SpkSortoolbox,'Units','Normalized');
    else
        Scrsize = get(0,'ScreenSize');
        fwidth=1010;
        fheight=min(Scrsize(4),960);
        toolboxpos=[floor(fwidth/2) floor(fheight/2) fwidth fheight];
    end


    fheight=480;
    fwidth=300;
    BGcolor=[0.7 0.7 0.7];
    Detectparamsdialog=figure('Name','Detection Params','NumberTitle','off','Visible','on','Position',[toolboxpos(1),toolboxpos(2),fwidth,fheight],'Menubar','none','Color',[0.7 0.7 0.7],...
        'CloseRequestFcn',{@(source,event)quitparams_Callback(source,event,obj)});
    txtsize=150;
    editsize=150;

    if isequal(obj.Params.fileinfo.electrode,[32])
        txtpopupelconfig=1;
    elseif isequal(obj.Params.fileinfo.electrode,[32 32])
        txtpopupelconfig=2;
    elseif isequal(obj.Params.fileinfo.electrode,[16])
        txtpopupelconfig=3;
    elseif isequal(obj.Params.fileinfo.electrode,[32 16])
        txtpopupelconfig=4;
    else
        txtpopupelconfig=5;
    end

    topDial=fheight-30;
    hTXTacqsystem = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','ACQ. SYSTEM',...
        'Position',[0.01,topDial,fwidth,20],'BackgroundColor',BGcolor);
    topDial=topDial-30;
    htxtacqsystem = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','aqcuisition system',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hacqsystem = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',obj.Params.fileinfo.acqsystem,...
        'Position',[txtsize,topDial,txtsize,20],'Callback',{@(source,event)acqsystem_Callback(source,event,obj)});

    topDial=topDial-30;
    htxtsamplingrate = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','Sampling rate',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hsamplingrate = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.fileinfo.samplingrate),...
        'Position',[txtsize,topDial,txtsize,20],'Callback',{@(source,event)samplingrate_Callback(source,event,obj)});

    topDial=topDial-30;
    hTXTelconfig = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','ELECTRODE',...
        'Position',[0.01,topDial,fwidth,20],'BackgroundColor',BGcolor);
    topDial=topDial-30;
    htxtChannelnum = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','Channels ID#',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hChannelnum = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.fileinfo.Channelnum),...
        'Position',[txtsize,topDial,txtsize,20],'Callback',{@(source,event)Channelnum_Callback(source,event,obj)});

    topDial=topDial-30;
    htxtelconfig = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','electrode configuration',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    helconfig = uicontrol('Parent',Detectparamsdialog,'Style','popupmenu',...
        'String',{'1 x 32 channels','2 x 32 channels','1 x 16 channels','32 + 16 channels','custom'},'Value',txtpopupelconfig,...
        'Position',[txtsize,topDial,txtsize,20],'Callback',{@(source,event)elecconfig_Callback(source,event,obj)});

    topDial=topDial-30;
    htxtFILTERING = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','FILTERING',...
        'Position',[0.01,topDial,fwidth,20],'BackgroundColor',BGcolor);
    topDial=topDial-30;
    htxtseglength = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','Segment size (sec)',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hseglength = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.detect.length2load),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)length2load_Callback(source,event,obj)});
    topDial=topDial-30;
    htxtlowcutoff = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','low cutoff',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hlowcutoff = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.detect.FcLow),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)lowcutoff_Callback(source,event,obj)});
    topDial=topDial-20;
    htxthighcutoff = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','high cutoff',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hhighcutoff = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.detect.FcHigh),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)highcutoff_Callback(source,event,obj)});

    topDial=topDial-30;
    htxtbutterorder = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','Filter order',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hbutterorder = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.detect.butterOrder),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)butterOrder_Callback(source,event,obj)});
    topDial=topDial-30;
    hFsubmedian = uicontrol('Parent',Detectparamsdialog,'Style','checkbox',...
        'String','median filtering','Min',0,'Max',1,'Value',obj.Params.detect.Fsubmedian,...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)Fsubmedian_Callback(source,event,obj)});

    topDial=topDial-30;
    htxtDETECTION = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','DETECTION',...
        'Position',[0.01,topDial,fwidth,20],'BackgroundColor',BGcolor);
    topDial=topDial-20;
    htxtzthresh = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','z-threshold',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hzthresh = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.detect.zthresh),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)zthresh_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtwintime = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','waveform length (ms)',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hwintime = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(2*obj.Params.detect.wtime2),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)wintime_Callback(source,event,obj)});

    topDial=topDial-20;
    htxtTcensoredfactor = uicontrol('Parent',Detectparamsdialog,'Style','text',...
        'String','Censored factor',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hTcensoredfactor = uicontrol('Parent',Detectparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.detect.Tcensoredfactor),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)Tcensoredfactor_Callback(source,event,obj)});

    topDial=topDial-20;
    hok = uicontrol('Parent',Detectparamsdialog,'Style','pushbutton',...
        'String','0K',...
        'Position',[0.01,topDial,txtsize,20],'Callback',{@(source,event)ok_Callback(source,event,obj)});
    hcancel = uicontrol('Parent',Detectparamsdialog,'Style','pushbutton',...
        'String','Cancel',...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)Cancel_Callback(source,event,obj)});

        function acqsystem_Callback(source,event,obj)
            obj.Params.fileinfo.acqsystem=get(source, 'String');
        end

        function samplingrate_Callback(source,event,obj)
            obj.Params.fileinfo.samplingrate=str2double(get(source,'String'));
        end

        function Channelnum_Callback(source,event,obj)
            obj.Params.fileinfo.Channelnum=str2num(get(source,'String'));
            obj.Params.fileinfo.Channelcount=numel(obj.Params.fileinfo.Channelnum);
        end

        function elecconfig_Callback(source,event,obj)
            str = get(source, 'String');
            val = get(source,'Value');
            switch val
                case 1
                    obj.Params.fileinfo.electrode=32;
                case 2
                    obj.Params.fileinfo.electrode=[32 32];
                case 3
                    obj.Params.fileinfo.electrode=[16];
                case 4
                    obj.Params.fileinfo.electrode=[32 16];
                case 5
                    option=[];
                    option.WindowStyle='modal';
                    option.Resize='on';
                    option.Interpreter='none';
                    elconf=inputdlg({'# channels (per electrode)'},'electrode configuration',1,{num2str(obj.Params.fileinfo.electrode)},option);
                    if ~isempty(char(elconf))
                        obj.Params.fileinfo.electrode=round(str2num(elconf{1}));
                    end
            end
        end

        function length2load_Callback(source,event,obj)
            obj.Params.detect.length2load=str2double(get(source,'String'));
        end

        function lowcutoff_Callback(source,event,obj)
            obj.Params.detect.FcLow=str2double(get(source,'String'));
        end

        function highcutoff_Callback(source,event,obj)
            obj.Params.detect.FcHigh=str2double(get(source,'String'));
        end

        function butterOrder_Callback(source,event,obj)
            obj.Params.detect.butterOrder=str2double(get(source,'String'));
        end

        function Fsubmedian_Callback(source,event,obj)
            val = get(source,'Value');
            if val==1
                obj.Params.detect.Fsubmedian=true;
            else
                obj.Params.detect.Fsubmedian=false;
            end
        end

        function zthresh_Callback(source,event,obj)
            obj.Params.detect.zthresh=str2double(get(source,'String'));
        end

        function wintime_Callback(source,event,obj)
            obj.Params.detect.wtime2=round(0.5*str2double(get(source,'String')));
            obj.Params.detect.wtime1=-obj.Params.detect.wtime2;
        end

        function Tcensoredfactor_Callback(source,event,obj)
            obj.Params.detect.Tcensoredfactor=str2double(get(source,'String'));
        end

        function quitparams_Callback(source,event,obj)
            delete(Detectparamsdialog);
            if ~ischar(obj.SpkSortoolbox)
                obj.UpdateParamsPanel;
            end
        end

        function ok_Callback(source,event,obj)
            delete(Detectparamsdialog);
            if ~ischar(obj.SpkSortoolbox)
                obj.UpdateParamsPanel;
            end
            obj.Spkdetection;
            obj.SaveSpikeData;
        end

        function Cancel_Callback(source,event,obj)
            delete(Detectparamsdialog);
            if ~ischar(obj.SpkSortoolbox)
                obj.UpdateParamsPanel;
            end
        end

    set([Detectparamsdialog,...
        hTXTacqsystem,htxtacqsystem,hacqsystem,htxtsamplingrate,hsamplingrate,...
        hTXTelconfig,htxtChannelnum,hChannelnum,htxtelconfig,helconfig,htxtseglength,hseglength...
        htxtFILTERING,htxtlowcutoff,hlowcutoff,htxthighcutoff,hhighcutoff,htxtbutterorder,hbutterorder,hFsubmedian,...
        htxtDETECTION,htxtzthresh,hzthresh,htxtwintime,hwintime,htxtTcensoredfactor,hTcensoredfactor,hok,hcancel],'Units','normalized');
end