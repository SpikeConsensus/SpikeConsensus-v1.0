function obj=InstallDialogcluster(obj)
    parentpanel=obj.dialogcluster;

    set(parentpanel,'Units','Pixels');
    wsize=get(parentpanel,'Position');
    set(parentpanel,'Units','normalized');
    dcolor=get(parentpanel,'Backgroundcolor');
    topDial=wsize(4)-0.05*wsize(4);
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

    Funit=false;
    Fdispoversampled=false;
    Fdisprealigned=false;
    Fdispwhitened=true;
    Fdispclass=true;
    Fdispcluttest=false;
    Fdispclusttree=false;
    Fdispstd=false;
    Fdispchauto=true;

    hFdispoversampled = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','oversample','Min',0,'Max',1,'Value',Fdispoversampled,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdispoversampled_Callback(source,event,obj)});

    topDial=topDial-30;
    hFdisprealigned = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','realign','Min',0,'Max',1,'Value',Fdisprealigned,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdisprealigned_Callback(source,event,obj)});

    topDial=topDial-30;
    hFdispmergedunit = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','Overlay','Min',0,'Max',1,'Value',Funit,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdispunit_Callback(source,event,obj)});              

    topDial=topDial-30;
    hFdispwhitened = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','Whitened','Min',0,'Max',1,'Value',Fdispwhitened,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdispwhitened_Callback(source,event,obj)});

    topDial=topDial-30;                          
    hFdispstd = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','disp. STD','Min',0,'Max',1,'Value',Fdispstd,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdispstd_Callback(source,event,obj)});
    topDial=topDial-30;                          
    hFdisptest = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','Post. Test','Min',0,'Max',1,'Value',Fdispcluttest,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdispcluttest_Callback(source,event,obj)});

    topDial=topDial-30;              
    htxtlistchannels = uicontrol('Parent',parentpanel,'Style','text',...
                                 'String','channels','Units','Pixels','Position',[rightDial,topDial,dwidth,20],'BackgroundColor',dcolor);

    topDial=topDial-150;                     
    hlistchannels=uicontrol('Parent',parentpanel,'Style','listbox','max',selmax,'min',selmin,...
                            'String',strch,'Units','Pixels','Position',[rightDial,topDial,dwidth,140],'BackgroundColor',dcolor,...
                            'Callback',{@(source,event)dispchcluster_Callback(source,event,obj)});                                
    topDial=topDial-40;                          
    hFdispchauto = uicontrol('Parent',parentpanel,'Style','checkbox',...
                      'String','auto ch. select','Min',0,'Max',1,'Value',Fdispchauto,...
                      'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Fdispchauto_Callback(source,event,obj)});

    topDial=topDial-30;                     
    htxtlistunits = uicontrol('Parent',parentpanel,'Style','text',...
                                 'String','Clusters','Units','Pixels','Position',[rightDial,topDial,dwidth,20],'BackgroundColor',dcolor);

    topDial=topDial-150;                         
    hlistunits=uicontrol('Parent',parentpanel,'Style','listbox','max',selmax,'min',selmin,...
                         'String',strunit,'Units','Pixels','Position',[rightDial,topDial,dwidth,140],'BackgroundColor',dcolor,...
                         'Callback',{@(source,event)dispunitcluster_Callback(source,event,obj)});


    topDial=topDial-30;                         
    hupdatePmiss=uicontrol('Parent',parentpanel,'Style','pushbutton','String','Update Pmiss',...
                         'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)UpdatePmiss_Callback(source,event,obj)});
    topDial=topDial-30;                         
    hshowdendro=uicontrol('Parent',parentpanel,'Style','checkbox',...
                          'String','Tree','Min',0,'Max',1,'Value',Fdispclusttree,...
                          'Position',[rightDial,topDial,dwidth,20],'Callback',{@(source,event)Showdendro_Callback(source,event,obj)});

    set([parentpanel,...
         hFdisprealigned,hFdispoversampled,hFdispmergedunit,hFdispwhitened,...
         hFdisptest,hFdispchauto,hFdispstd,htxtlistchannels,hlistchannels,...
         htxtlistunits,hlistunits,hupdatePmiss,hshowdendro],'Units','normalized');

    function Fdispoversampled_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispoversampled=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        else
            Fdispoversampled=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        end
    end

    function Fdisprealigned_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdisprealigned=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        else
            Fdisprealigned=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        end
    end

    function Fdispunit_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Funit=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        else
            Funit=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        end
    end

    function Fdispwhitened_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispwhitened=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        else
            Fdispwhitened=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        end
    end

    function Fdispclass_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispclass=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);                    
        else
            Fdispclass=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);
        end
    end

    function Fdispcluttest_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispcluttest=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);                   
        else
            Fdispcluttest=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);
        end
    end

    function Fdispchauto_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispchauto=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);                  
        else
            Fdispchauto=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);
        end
    end                        

    function Fdispstd_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispstd=true;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);
        else
            Fdispstd=false;
            obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
            obj.clustprofile(UnitSel);
        end
    end            

    function dispchcluster_Callback(source,event,obj)                
        str=get(source,'string');
        val=get(source,'Value');
        ChSel=round(str2double(str(val)))';
        obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);                
        obj.clustprofile(UnitSel);
    end

    function dispunitcluster_Callback(source,event,obj)
        str=get(source,'string');
        val=get(source,'Value');
        UnitSel=round(str2double(str(val)))';
        clustSeldisp = [];
        cSelAll = [];
        for c = 1:numel(UnitSel)
            uSel=unique(obj.getUnitIDMultiCh(UnitSel(c)));
            cSel=[];
            for u=1:numel(uSel)
                cSel=[cSel obj.getClusterIDMultiCh(uSel(u))];
            end
            cSelAll = [cSelAll cSel];            
            clustSel = sort(unique(cSel),'ascend');
            idx = find(clustSel == UnitSel(c));
            clustSeldisp = [clustSeldisp clustSel(max(1, idx-2):min(idx+2, end))];
        end
        obj.editUnitID(unique(cSelAll)');
        Updateclustertree(obj, unique(cSelAll)');
        UnitSel = clustSeldisp;
        obj.ShowClusterquality(unique(UnitSel));
        obj.ShowUnitquality(unique(UnitSel));
        obj.Updatepagecluster(ChSel,UnitSel,Funit,Fdispoversampled,Fdisprealigned,Fdispwhitened,Fdispclass,Fdispstd,Fdispchauto);
        obj.clustprofile(UnitSel);
    end 

    function UpdatePmiss_Callback(source,event,obj)
        UpdateTotalMiss(obj,true);
    end

    function Showdendro_Callback(source,event,obj)
        val = get(source,'Value');
        if val==1
            Fdispclusttree=true;
            Updateclustertree(obj,UnitSel);
        else
            Fdispclusttree=false;
            obj.UpdateParamsPanel;
        end
    end

    function Updateclustertree(obj, UnitSel)
        el=floor((obj.MclusterChID(UnitSel(1))-1)/obj.Params.fileinfo.electrode(1))+1;
        clust2plot=find(ismember(obj.MclusterChID,obj.elconfig(1,el):obj.elconfig(2,el)));
        if isfield(obj.Quality,'ConfusionCMat')
            rng(1);
            distUmat=squareform(1-obj.Quality.ConfusionCMat(clust2plot,clust2plot));
            nbclust=size(obj.Quality.ConfusionCMat(clust2plot,clust2plot),1);
            Z=linkage(distUmat);
            f=figure;
            try
                %                     leafOrder=optimalleaforder(Z,distUmat,'Criteria','group');
                %                     [~,~,outperm]=dendrogram(Z,nbclust,'Reorder',leafOrder);%dendrogram(Z,nbclust,'Reorder',[1:23]);%dendrogram(Z,nbclust,'Reorder',leafOrder);
                [H,~,outperm]=dendrogram(Z,nbclust);
            catch
                [H,~,outperm]=dendrogram(Z,nbclust);
            end
            if ishandle(obj.paramspanel)
                delete(get(obj.paramspanel,'Children'));
            end
            axs=axes('Parent',obj.paramspanel);
            set(H,'Parent',axs);
            set(axs,'XTick',1:nbclust);
            set(axs,'XTickLabel',clust2plot(outperm)');
            selidx = find(ismember(clust2plot(outperm), UnitSel));
            hold(axs, 'on');
            plot(axs, selidx, ones(size(selidx)) + 0.02, '+r');
            hold(axs, 'off');
            set(axs,'Ylim',[0.7 1.05],'TickLength',[0 0],'FontSize',5);
            delete(f);            
        end
    end
end