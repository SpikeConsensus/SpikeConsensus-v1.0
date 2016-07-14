function obj=SpikeSortToolbox(obj)
    Scrsize = get(0,'ScreenSize');
    obj.SpkSortoolbox=figure('Name','SpikeSorting Toolbox','NumberTitle','off','Visible','on','Position',[0,0,1280,min(Scrsize(4),960)],'Menubar','none','Color',[1 1 1],...
        'CloseRequestFcn',{@(source,event)quit_Callback(source,event,obj)});
    obj.loadmenu = uimenu('Parent',obj.SpkSortoolbox,'Label','File');
    uimenu('Parent',obj.loadmenu,'Label','load file','Callback',{@(source,event)LoadData_Callback(obj,source,event)});
    uimenu('Parent',obj.loadmenu,'Label','save file','Callback',{@(source,event)SaveData_Callback(obj,source,event)});
    uimenu('Parent',obj.loadmenu,'Label','Append file','Callback',{@(source,event)Appendfile_Callback(obj,source,event)});
    uimenu('Parent',obj.loadmenu,'Label','undo Append file','Callback',{@(source,event)resetAppendfile_Callback(obj,source,event)});

    obj.detectionmenu = uimenu('Parent',obj.SpkSortoolbox,'Label','Detection');
    uimenu('Parent',obj.detectionmenu,'Label','execute','Callback',{@(source,event)dodetection_Callback(obj,source,event)});

    obj.clustermenu = uimenu('Parent',obj.SpkSortoolbox,'Label','Sorting');
    uimenu('Parent',obj.clustermenu,'Label','Load clusters','Callback',{@(source,event)LoadSortedData_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Clustering params','Callback',{@(source,event)MClusterparams_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Test NbKcluster space','Callback',{@(source,event)NbKclustertest_Callback(obj,source,event)});    
    uimenu('Parent',obj.clustermenu,'Label','Clustering (iteration phase)','Callback',{@(source,event)iteCluster_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Clustering (consensus phase)','Callback',{@(source,event)ConsensusCluster_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Clustering (classification phase)','Callback',{@(source,event)ClassifyCluster_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Quality measures','Callback',{@(source,event)QualityCluster_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Reorder clusters','Callback',{@(source,event)ReorderCluster_Callback(obj,source,event)});
    uimenu('Parent',obj.clustermenu,'Label','Save clusters','Callback',{@(source,event)SaveSortedData_Callback(obj,source,event)});

    fwidth=1;
    fheight=1;
    obj.dialogClustquality = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.005*fwidth,0.06*fheight,0.15*fwidth,0.60*fheight],'BackgroundColor',[0.8 0.8 0.8]);
    obj.dialogunit = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.005*fwidth,0.67*fheight,0.15*fwidth,0.20*fheight],'BackgroundColor',[0.8 0.8 0.8]);
    obj.dialogprofile = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.85*fwidth,0.67*fheight,0.15*fwidth,0.20*fheight],'BackgroundColor',[0.8 0.8 0.8]);
    obj.dialogUnitquality = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.85*fwidth,0.06*fheight,0.15*fwidth,0.60*fheight],'BackgroundColor',[0.8 0.8 0.8]);

    obj.dialogdetect = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.15*fwidth,0.06*fheight,0.05*fwidth,0.81*fheight],'BackgroundColor',[1 1 1],'Visible','off');
    obj.dialogcluster = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.15*fwidth,0.06*fheight,0.05*fwidth,0.81*fheight],'BackgroundColor',[1 1 1],'Visible','off');
    obj.dialogclassify = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.15*fwidth,0.06*fheight,0.05*fwidth,0.81*fheight],'BackgroundColor',[1 1 1],'Visible','off');

    obj.pagedetect = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.2*fwidth,0.06*fheight,0.65*fwidth,0.81*fheight],'BackgroundColor',[1 1 1],'Visible','off');
    obj.pagecluster = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.2*fwidth,0.06*fheight,0.65*fwidth,0.81*fheight],'BackgroundColor',[1 1 1],'Visible','off');
    obj.pageclassify = uipanel('Parent',obj.SpkSortoolbox,'Units','normalized','Position',[0.2*fwidth,0.06*fheight,0.65*fwidth,0.81*fheight],'BackgroundColor',[1 1 1],'Visible','off');

    obj.pagebuttondetect = uicontrol('Parent',obj.SpkSortoolbox,'Style','togglebutton','String','DETECTION','Units','normalized',...
        'Position',[0.005*fwidth,0.005*fheight,0.33*fwidth,0.05*fheight],...
        'Callback',{@(source,event)pagedetect_Callback(source,event,obj)});
    obj.pagebuttoncluster = uicontrol('Parent',obj.SpkSortoolbox,'Style','togglebutton','String','CLUSTERS','Units','normalized',...
        'Position',[0.33*fwidth,0.005*fheight,0.33*fwidth,0.05*fheight],...
        'Callback',{@(source,event)pagecluster_Callback(source,event,obj)});
    obj.pagebuttonclassify = uicontrol('Parent',obj.SpkSortoolbox,'Style','togglebutton','String','CLASSIFICATION','Units','normalized',...
        'Position',[0.66*fwidth,0.005*fheight,0.33*fwidth,0.05*fheight],...
        'Callback',{@(source,event)pageclassify_Callback(source,event,obj)});
    obj=UpdateParamsPanel(obj);
    obj=InstallDialogdetect(obj);
    obj=InstallDialogcluster(obj);
    obj=InstallDialogclassify(obj);

    set([obj.SpkSortoolbox,...
        obj.pagedetect,obj.pagecluster,obj.pageclassify,...
        obj.dialogdetect,obj.dialogcluster,obj.dialogclassify,...
        obj.dialogunit,obj.dialogClustquality,obj.dialogprofile,obj.dialogUnitquality,...
        obj.pagebuttondetect,obj.pagebuttoncluster,obj.pagebuttonclassify],'Units','Normalized');
    
    function quit_Callback(source,event,obj)
        selection = questdlg('Are you sure you want to quit?',...
                'Close Request Function',...
                'Yes','No','Yes');
        if strcmp(selection,'Yes')
            evalin('base','clear SpikeSorter;');
            delete(obj.SpkSortoolbox);
        end
    end 
end

function pagedetect_Callback(source,event,obj)
    set([obj.pagecluster obj.dialogcluster obj.pageclassify obj.dialogclassify],'Visible','off');
    set([obj.pagebuttoncluster obj.pagebuttonclassify],'Value',0);
    set([obj.pagedetect obj.dialogdetect],'Visible','on');
    set(obj.pagebuttondetect,'Value',1);
end

function pagecluster_Callback(source,event,obj)
    set([obj.pagedetect obj.dialogdetect obj.pageclassify obj.dialogclassify],'Visible','off');
    set([obj.pagebuttondetect obj.pagebuttonclassify],'Value',0);
    set([obj.pagecluster obj.dialogcluster],'Visible','on');
    set(obj.pagebuttoncluster,'Value',1);
end

function pageclassify_Callback(source,event,obj)
    set([obj.pagedetect obj.dialogdetect obj.pagecluster obj.dialogcluster],'Visible','off');
    set([obj.pagebuttondetect obj.pagebuttoncluster],'Value',0);
    set([obj.pageclassify obj.dialogclassify],'Visible','on');
    set(obj.pagebuttonclassify,'Value',1);
end

function dodetection_Callback(obj,source,event)
    obj.executedetection(true);
end

function MClusterparams_Callback(obj,source,event)
    obj.SpkClusterdialog();
end                

function iteCluster_Callback(obj,source,event)
    obj.SpkClusterMultiCh();           
end

function NbKclustertest_Callback(obj,source,event)
    if exist([obj.datafolder filesep obj.filename '_chi2vsKite.mat']);
        choice = questdlg('already done. Wanna do it again?','chi2vsKite','Yes','No','No');
    else
        choice = 'Yes';
    end
    if strcmp(choice,'Yes')
        obj.SpkClusterMultiCh(true);
    else
        load([obj.datafolder filesep obj.filename '_chi2vsKite.mat'],'chi2_95','nbunitsglobal0');
        figure;
        plot(nbunitsglobal0,chi2_95);
        hold on;
        plot(nbunitsglobal0(2:end),diff(chi2_95),'r');
    end
end

function ConsensusCluster_Callback(obj,source,event)
    obj.SpkMetaClusterMultiCh();
end

function ClassifyCluster_Callback(obj,source,event)
    obj.SpkMetaClassifyMultiCh();
end

function QualityCluster_Callback(obj,source,event)
    obj.SpkMetaUnitQuality();
    obj.SpkMetaUnitQuality(true);
end

function ReorderCluster_Callback(obj,source,event)
    obj.ReorderClusters();
end

function SaveSortedData_Callback(obj,source,event)
    obj.SaveSortedData;
end                        

function LoadSortedData_Callback(obj,source,event)
    obj.LoadSortedSpikes('',true);
end

function LoadData_Callback(obj,source,event)
    obj.LoadFileParams();
    obj.LoadSpikeData();
    obj.LoadSortedSpikes();
    obj.UpdateParamsPanel;
end

function SaveData_Callback(obj,source,event)
    obj.SaveSortedSpikes();
end

function Appendfile_Callback(obj,source,event)
    obj.Appendfile();
end

function resetAppendfile_Callback(obj,source,event)
    obj.resetAppendfile();
end