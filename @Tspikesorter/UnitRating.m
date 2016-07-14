function obj=UnitRating(obj)
    nbunit=obj.nbUnitMultiCh;
    nbUCol = 70;
    figrating=figure('Name','Unit Rating','NumberTitle','off','Visible','on','units','normalized','Position',[0,0,0.15+0.15*(nbunit>nbUCol),1],'Menubar','none','Color',[0.7 0.7 0.7]);
    BGcolor=[0.8 0.8 0.8];    
    if ~isfield(obj.Quality,'Rating') || size(obj.Quality.Rating,1) ~= nbunit
        obj.Quality.Rating=false(nbunit,4);
    end
    clustunitID=obj.Quality.Rating;

    txtsize=1/(6+6*(nbunit>nbUCol));
    txtheight=0.05;
    editsize=0.8;
    topDial=0.95;
    htxt0= uicontrol('Parent',figrating,'Style','text',...
                     'String','# Pmis SNR ISI All','Units','Normalized',...
                     'Position',[0.01,topDial,1,txtheight],'BackgroundColor',BGcolor);
    hedit=cell(nbunit,4);
    for u=1:nbunit
        for i=1:4
            htxt1= uicontrol('Parent',figrating,'Style','text',...
                             'String',num2str(u),'Units','Normalized',...
                             'Position',[0.01 + 0.5*(u>nbUCol),topDial-(u - nbUCol*(u>nbUCol))*1/(min(nbUCol,nbunit)+5),txtsize,1/(min(nbUCol,nbunit)+5)],'BackgroundColor',BGcolor);
            hedit{u,i} = uicontrol('Parent',figrating,'Style','checkbox',...
                                   'String','','Min',0,'Max',1,'Value',clustunitID(u,i),'Units','Normalized',...
                                   'Position',[0.01 + 0.5*(u>nbUCol)+txtsize*i,topDial-(u - nbUCol*(u>nbUCol))*1/(min(nbUCol,nbunit)+5),txtsize,1/(min(nbUCol,nbunit)+5)],'Callback',{@(source,event)unitRating_Callback(source,event,u,i,obj)});                                  
        end
    end            

    function unitRating_Callback(source,event,u,i,obj)
        val = get(source,'Value');
        obj.Quality.Rating(u,i)=val;
    end                                
end