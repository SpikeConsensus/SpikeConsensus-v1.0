function obj=editUnitID(obj,clustSel)
    delete(get(obj.dialogunit,'Children'));
    BGcolor=[0.8 0.8 0.8];            

    clustunitID=obj.SpkclustIDMulti;

    txtsize=0.15;
    txtheight=0.08;
    editsize=0.5;
    topDial=0.92;
    htxt0= uicontrol('Parent',obj.dialogunit,'Style','text',...
                     'String','cluster#','Units','Normalized',...
                     'Position',[0.01,topDial,txtsize,txtheight],'BackgroundColor',BGcolor);
    htxt=cell(1,size(clustunitID,2));
    hedit=cell(1,size(clustunitID,2));
    hpmissedit=cell(1,size(clustunitID,2));
    topDial = 0.90;
    for u=1:size(clustSel,1)                
        htxt{u}= uicontrol('Parent',obj.dialogunit,'Style','text',...
                           'String',num2str(clustSel(u)),'Units','Normalized',...
                           'Position',[0.01,topDial-u*1/15,txtsize,txtheight],'BackgroundColor',BGcolor);
    end
    for u=1:size(clustSel,1) 
        hedit{u} = uicontrol('Parent',obj.dialogunit,'Style','edit',...
                             'String',num2str(clustunitID(clustSel(u),clustunitID(clustSel(u),:)>0)),'Units','Normalized',...
                             'Position',[txtsize,topDial-u*1/15,editsize,txtheight],'Callback',{@(source,event)unitID_Callback(source,event,clustSel(u),obj)});
        hpmissedit{u} = uicontrol('Parent',obj.dialogunit,'Style','edit',...
                                  'String',num2str(obj.SpkclustPmisMulti(clustSel(u))),'Units','Normalized',...
                                  'Position',[txtsize + editsize,topDial-u*1/15,0.2,txtheight],'Callback',{@(source,event)pmisunit_Callback(source,event,u,obj)});    
    end            

    function unitID_Callback(source,event,u,obj)
        ID = str2num(get(source,'String')); 
        obj.SpkclustIDMulti(u,:) = zeros(1,size(obj.SpkclustIDMulti,2));
        if size(ID,2)>0
            obj.SpkclustIDMulti(u,1:size(ID,2)) = ID;
        end                
    end 
    function pmisunit_Callback(source,event,u,obj)
        pmis = str2num(get(source,'String'));
        
        distUmat=squareform((1-obj.Quality.ConfusionCMat));
        Z=linkage(distUmat);
        metaClustID=cluster(Z,'cutoff',1-pmis,'criterion','distance');
        
        clustSelID = metaClustID(clustSel(u));
        IDs = find(metaClustID == clustSelID);
        IDs = IDs(:)';
        
%         IDs = find(obj.Quality.ConfusionCMat(clustSel(u),:) > pmis);
        [r, ~] = find(ismember(obj.SpkclustIDMulti, IDs));
        r(r == clustSel(u)) = [];
        IDtemp = obj.SpkclustIDMulti(r,:);
        IDtemp = IDtemp(IDtemp > 0);
        IDtemp = IDtemp(:)';
        IDs = unique([IDs unique(IDtemp)]);
        
        obj.SpkclustPmisMulti(IDs) = pmis;%instead outsource to another function which adjust all Pmiss of all concerned clusters
        restore = obj.SpkclustIDMulti(IDs(IDs == clustSel(u)),:);
        obj.SpkclustIDMulti(IDs(IDs == clustSel(u)),1:numel(IDs)) = IDs;
        obj.SpkclustIDMulti(IDs(IDs ~= clustSel(u)),:) = 0;        
        obj.SpkclustIDMulti(IDs(IDs == clustSel(u)),numel(IDs)+1:end) = 0;        
        restore = restore(restore > 0 & ~ismember(restore,IDs));
        obj.SpkclustIDMulti(restore,1) = restore;
        obj.SpkclustIDMulti(restore,2:end) = 0;
        obj.SpkclustPmisMulti(restore) = pmis;
        
        for clust = 1:numel(restore)
            clustSelID = metaClustID(restore(clust));
            IDs2 = find(metaClustID == clustSelID);
            IDs2 = IDs2(:)';
            
            %         IDs = find(obj.Quality.ConfusionCMat(clustSel(u),:) > pmis);
            [r, ~] = find(ismember(obj.SpkclustIDMulti, IDs2));
            r(r == restore(clust)) = [];
            IDtemp = obj.SpkclustIDMulti(r,:);
            IDtemp = IDtemp(IDtemp > 0);
            IDtemp = IDtemp(:)';
            IDs2 = unique([IDs2 unique(IDtemp)]);
            
            obj.SpkclustPmisMulti(IDs2) = pmis;%instead outsource to another function which adjust all Pmiss of all concerned clusters
%             restore = obj.SpkclustIDMulti(IDs(IDs == restore(clust)),:);
            obj.SpkclustIDMulti(IDs2(IDs2 == restore(clust)),1:numel(IDs2)) = IDs2;
            obj.SpkclustIDMulti(IDs2(IDs2 ~= restore(clust)),:) = 0;
%             obj.SpkclustIDMulti(IDs(IDs == restore(clust)),numel(IDs)+1:end) = 0;
%             restore = restore(restore > 0 & ~ismember(restore,IDs));
%             obj.SpkclustIDMulti(restore,1) = restore;
%             obj.SpkclustIDMulti(restore,2:end) = 0;
%             obj.SpkclustPmisMulti(restore) = pmis;
        end
        
        for u=1:size(clustSel,1)
            if ismember(clustSel(u),[IDs restore])
                set(hpmissedit{u}, 'String', num2str(pmis));
            end
            set(hedit{u}, 'String', num2str(obj.SpkclustIDMulti(clustSel(u),obj.SpkclustIDMulti(clustSel(u),:)>0)));
        end
    end 
end