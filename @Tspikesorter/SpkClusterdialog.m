function obj = SpkClusterdialog(obj)
    if ~ischar(obj.SpkSortoolbox)
        set(obj.SpkSortoolbox,'Units','Pixels');
        toolboxpos=get(obj.SpkSortoolbox,'Position');
        set(obj.SpkSortoolbox,'Units','Normalized');
    else
        Scrsize = get(0,'ScreenSize');
        fwidth=1010;
        fheight=min(Scrsize(4),960);
        toolboxpos=[0 0 fwidth fheight];
    end

    fheight=380;
    fwidth=300;
    BGcolor=[0.7 0.7 0.7];
    Clusterparamsdialog=figure('Name','Clustering Params','NumberTitle','off','Visible','on','Position',[toolboxpos(1),toolboxpos(2),fwidth,fheight],'Menubar','none','Color',[0.7 0.7 0.7],...
        'CloseRequestFcn',{@(source,event)quitparams_Callback(source,event,obj)});
    txtsize=150;
    editsize=150;
    topDial=fheight-0.1*fheight;
    htxtCLUSTERING = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','CLUSTERING & CLASSIFICATION',...
        'Position',[0.01,topDial,fwidth,20],'BackgroundColor',BGcolor);
    topDial=topDial-20;
    htxtPCAcompMax = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','max # PC per channel',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hPCAcompMax = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.nbPCmax),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)PCAcompMax_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtPCAcompMin = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','min # PC per channel',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hPCAcompMin = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.nbPCmin),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)PCAcompMin_Callback(source,event,obj)});
    topDial=topDial-20;
    htxttrainingsmaxnbspk = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','max # clean Spikes',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    htrainingsmaxnbspk = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.TrainingMaxnbSpk),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)TrainingMaxnbSpk_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtnbiteration = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','# of iterations',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hnbiteration = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.NbIteration),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)NbIteration_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtnbKcluster = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','# of Kclusters',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hnbKcluster = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.NbKcluster),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)NbKcluster_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtnbKmeansRep = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','# of Kmeans replicates',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hnbKmeansRep = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.KmeansReplicates),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)NbKmeansReplicates_Callback(source,event,obj)});
    topDial=topDial-20;
    hFKmeansOnline = uicontrol('Parent',Clusterparamsdialog,'Style','checkbox',...
        'String','Kmeans online updates','Min',0,'Max',1,'Value',strcmp(obj.Params.MCluster.FKmeansOnline, 'on'),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)FKmeansOnline_Callback(source,event,obj)});    
    topDial=topDial-20;
    htxtseed = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','initial Seed',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hseed = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.Seed),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)Seed_Callback(source,event,obj)});   
    topDial=topDial-20;
    htxtnbmaxcluster = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','max #clusters per ch.',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hnbmaxcluster = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.nbmaxunits),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)nbmaxcluster_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtprojrange = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Projection range',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hprojrange = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.ProjRange),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)ProjRange_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtCoreClustersMeth = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Core clusters Method',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    
    if strcmp(obj.Params.MCluster.CoreClustersMeth,'Consistency')
        val = 1;
    elseif strcmp(obj.Params.MCluster.CoreClustersMeth,'Best')
        val = 2;
    else
        val = 3;
    end
    hCoreClustersMeth = uicontrol('Parent',Clusterparamsdialog,'Style','popupmenu',...
        'String',{'Consistency','Best','not an option'},'Value',val,...
        'Position',[txtsize,topDial,txtsize,20],'Callback',{@(source,event)CoreClustersMeth_Callback(source,event,obj)});
    
    topDial=topDial-20;
    htxtminusize = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Cluster size thresh. values',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hminusize = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.minUsize),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)minUsize_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtPmisthreshold = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Pmis threshold',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hPmisthreshold = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.Pmisthreshold),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)Pmisthreshold_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtcdfChi2thresh = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Chi2 threshold (on cdf)',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hcdfChi2thresh = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.cdfChi2thresh),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)cdfChi2thresh_Callback(source,event,obj)});
    topDial=topDial-20;
    htxtcdfChi2threshnoise = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Chi2 threshold noise',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hcdfChi2threshnoise = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.cdfChi2threshnoise),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)cdfChi2threshnoise_Callback(source,event,obj)});    
    topDial=topDial-20;
    htxtUpdateTime = uicontrol('Parent',Clusterparamsdialog,'Style','text',...
        'String','Update fit (sec)',...
        'Position',[0.01,topDial,txtsize,20],'HorizontalAlignment','left','BackgroundColor',BGcolor);
    hUpdateTime = uicontrol('Parent',Clusterparamsdialog,'Style','edit',...
        'String',num2str(obj.Params.MCluster.UpdateTime),...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)UpdateTime_Callback(source,event,obj)});
    topDial=topDial-20;
    hFgreedy = uicontrol('Parent',Clusterparamsdialog,'Style','checkbox',...
        'String','Greedy method','Min',0,'Max',1,'Value',obj.Params.MCluster.Fgreedy,...
        'Position',[txtsize,topDial,editsize,20],'Callback',{@(source,event)Fgreedy_Callback(source,event,obj)});


        function PCAcompMax_Callback(source,event,obj)
            obj.Params.MCluster.nbPCmax=str2double(get(source,'String'));
        end
    
        function PCAcompMin_Callback(source,event,obj)
            obj.Params.MCluster.nbPCmin=str2double(get(source,'String'));
        end

        function nbmaxcluster_Callback(source,event,obj)
            obj.Params.MCluster.nbmaxunits=str2double(get(source,'String'));
        end

        function TrainingMaxnbSpk_Callback(source,event,obj)
            obj.Params.MCluster.TrainingMaxnbSpk=str2double(get(source,'String'));
        end

        function NbIteration_Callback(source,event,obj)
            obj.Params.MCluster.NbIteration=str2double(get(source,'String'));
        end
    
        function NbKcluster_Callback(source,event,obj)
            obj.Params.MCluster.NbKcluster=str2double(get(source,'String'));
        end
    
        function NbKmeansReplicates_Callback(source,event,obj)
            obj.Params.MCluster.KmeansReplicates=str2double(get(source,'String'));
        end
    
        function FKmeansOnline_Callback(source,event,obj)
            val = get(source,'Value');
            if val==1
                obj.Params.MCluster.FKmeansOnline='on';
            else
                obj.Params.MCluster.FKmeansOnline='off';
            end
        end            
    
        function Seed_Callback(source,event,obj)
            obj.Params.MCluster.Seed=str2double(get(source,'String'));
        end            

        function ProjRange_Callback(source,event,obj)
            obj.Params.MCluster.ProjRange=str2double(get(source,'String'));
        end
    
        function CoreClustersMeth_Callback(source,event,obj)
            str = get(source, 'String');
            val = get(source,'Value');
            switch val
                case 1
                    obj.Params.MCluster.CoreClustersMeth = 'Consistency';
                case 2
                    obj.Params.MCluster.CoreClustersMeth = 'Best';                
            end
        end            

        function minUsize_Callback(source,event,obj)
            obj.Params.MCluster.minUsize=str2num(get(source,'String'));
        end

        function Pmisthreshold_Callback(source,event,obj)
            obj.Params.MCluster.Pmisthreshold=str2double(get(source,'String'));
        end

        function cdfChi2thresh_Callback(source,event,obj)
            obj.Params.MCluster.cdfChi2thresh=str2double(get(source,'String'));
        end
    
        function cdfChi2threshnoise_Callback(source,event,obj)
            obj.Params.MCluster.cdfChi2threshnoise=str2double(get(source,'String'));
        end

        function UpdateTime_Callback(source,event,obj)
            obj.Params.MCluster.UpdateTime=str2double(get(source,'String'));
        end

        function Fgreedy_Callback(source,event,obj)
            val = get(source,'Value');
            if val==1
                obj.Params.MCluster.Fgreedy=true;
            else
                obj.Params.MCluster.Fgreedy=false;
            end
        end

        function quitparams_Callback(source,event,obj)
            delete(Clusterparamsdialog);
            if ~ischar(obj.SpkSortoolbox)
                obj.UpdateParamsPanel;
            end
        end

    set([Clusterparamsdialog,htxtCLUSTERING,...
        htxtPCAcompMax,hPCAcompMax,htxtPCAcompMin,hPCAcompMin,htxtnbmaxcluster,hnbmaxcluster,htxttrainingsmaxnbspk,htrainingsmaxnbspk,htxtnbiteration,hnbiteration,htxtnbKcluster,hnbKcluster,...
        htxtnbKmeansRep,hnbKmeansRep,hFKmeansOnline,htxtseed,hseed,htxtprojrange,hprojrange,htxtCoreClustersMeth,hCoreClustersMeth,htxtminusize,hminusize,htxtPmisthreshold,hPmisthreshold,...
        htxtcdfChi2thresh,hcdfChi2thresh,htxtcdfChi2threshnoise,hcdfChi2threshnoise,htxtUpdateTime,hUpdateTime,hFgreedy],'Units','normalized');
end