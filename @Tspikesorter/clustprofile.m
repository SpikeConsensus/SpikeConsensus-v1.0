function obj=clustprofile(obj,Usel)
    delete(get(obj.dialogprofile,'Children'));
    clustprofilemat=cell(size(Usel,2),1);
    hax=zeros(size(Usel,2),1); 
    hplot=zeros(size(Usel,2),1); 
    for u=1:size(Usel,2)
        hax(u)=subplot(1,size(Usel,2),u,'Parent',obj.dialogprofile,'XTickMode','manual','XTickLabel',[]);                
    end
    zmax=0.0001;
    for u=1:size(Usel,2)
        if Usel(u)<=size(obj.MclusterChID,2)
            uClustAve=obj.SpkwaveclustAve{Usel(u)};
            wnbpts=size(uClustAve,2);
            ChPosX=obj.Params.fileinfo.ChannelPosX;%round(obj.Params.fileinfo.ChannelPosX / max(diff(obj.Params.fileinfo.ChannelPosX))) + 1; %-min(obj.Params.fileinfo.ChannelPosX)+1;
            ChPosY=obj.Params.fileinfo.ChannelPosY;%round(obj.Params.fileinfo.ChannelPosY / max(diff(obj.Params.fileinfo.ChannelPosY))) + 1; %-min(obj.Params.fileinfo.ChannelPosY)+1;
            ChPosX(ChPosX>max(ChPosY)*10)=2;            
            clustprofilemat{u}=zeros(max(ChPosY),max(ChPosX));%zeros(size(uClustAve,1)+1,2);
            for ch=1:size(uClustAve,1)
                clustprofilemat{u}(ChPosY(ch),ChPosX(ch))=(sum(uClustAve(ch,floor(wnbpts/2)+1).^2)).^0.5;
            end
            if max(max(abs(clustprofilemat{u})))>zmax
                zmax=max(max(abs(clustprofilemat{u})));
            end
            hplot(u)=imagesc(clustprofilemat{u},'Parent',hax(u));
            set(get(hax(u),'Title'),'String',['#' num2str(Usel(u))]);
            if u>1
                set(hax(u),'YTickMode','manual','YTickLabel',[]);
            end   
        end
    end
    for u=1:size(Usel,2)
        set(hax(u),'CLim',[0 zmax]);
    end            
end