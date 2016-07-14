function obj=Updatepagecluster(obj,ChSel,UnitSel,Funit,Foversampled,Frealign,Fwhitened,Fsorted,Fdispstd,Fdispchauto)
    parentpanel=obj.pagecluster;   

    Fsorted = true;

    if Fwhitened
        maxY=30;
    else
        maxY=1000;%150;%150;%150;%150;%40;%150;
    end

%             Wn=[300 3000]/(obj.Params.fileinfo.samplingrate/2);
%             [butHB,butHA]=butter(6,Wn(1),'high');
%             [butLB,butLA]=butter(6,Wn(2),'low');
%             
    dispUnitSel=cell(1,size(UnitSel,2));
    if ~Funit
        for u=1:size(UnitSel,2)
            dispUnitSel{u}=UnitSel(u);
        end                
    else
        for u=1:size(UnitSel,2)
%                     [chref bla]=find(clustID==UnitSel(u));
%                     [unitch bla]=find(squeeze(clustID(chref,:,:))==UnitSel(u));
%                     dispUnitSel{u}=squeeze(clustID(chref,unitch,clustID(chref,unitch,:)>0))';
            dispUnitSel{u}=UnitSel;
        end                
    end

    if ~Fsorted
        clustens=obj.SpkwaveclustMulti;
    else
        %clustens=obj.SpkwaveclassMulti;
        unitdisp=unique(cell2mat(dispUnitSel));
        clustens=cell(max(unitdisp),1);
        %spkwave=obj.NoiseWhitening();
        %IDXglobalbisAll0=evalin('base','IDXglobalbisAll0');
        for u=1:size(unitdisp,2)
%                     if u==1
%                         spktimeU=find(ismember(obj.Spkevent(6,:),[-1]) & obj.Spkevent(10,:)<=2);
%                     else
            %spktimeU=find(ismember(IDXglobalbisAll0(:,1),unitdisp(u))' & obj.Spkevent(10,:)<=2);
            spktimeU=find(ismember(obj.Spkevent(6,:),unitdisp(u)) & obj.Spkevent(10,:)<=5);%& ((obj.Spkevent(2,:)-obj.Spkevent(2,1))*10^-3)>260*60000 & ((obj.Spkevent(2,:)-obj.Spkevent(2,1))*10^-3)<310*60000);
%                     end
            nbspk=min(300,numel(spktimeU));%numel(spktimeU);
            spktimeU=spktimeU(1:nbspk);

            clustens{unitdisp(u)}=obj.Spkwave(:,spktimeU(:,:),:);                    
        end
    end

    nbunit=size(UnitSel,2);

    chwinsize=3;
    chwin=cell(1,obj.Params.fileinfo.Channelcount);
    for el=1:size(obj.elconfig,2)
        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        nbCh=length(chfocus);
        for ch=1:nbCh
            chwin{chfocus(ch)}=find(abs(obj.Params.fileinfo.ChannelPosY-obj.Params.fileinfo.ChannelPosY(chfocus(ch)))<=chwinsize & abs(obj.Params.fileinfo.ChannelPosX-obj.Params.fileinfo.ChannelPosX(chfocus(ch)))<=chwinsize);    
        end
    end

    if Fdispchauto
        ChSel=[];
        for u=1:nbunit
            for c=1:size(dispUnitSel{u},2)
                ChSel=[ChSel chwin{round(obj.MclusterChID(dispUnitSel{u}(c)))}];
            end
        end
        ChSel=unique(ChSel(ChSel>=1 & ChSel<=size(obj.Spkwave,1)));
    end
%             ChSel=(ChSel*2-1)*2-1;
    nbCh=size(ChSel,2);            
    hwaveforms=zeros(nbCh,nbunit);            

    if ~isempty(get(parentpanel,'Children'))
        delete(get(parentpanel,'Children'));
    end

    for u=1:nbunit
        hwaveforms(1,u)=subplot(1,nbunit,u,'Parent',parentpanel);                    
    end                                                            

    if ~isempty(clustens{dispUnitSel{1}(1)})
        if Foversampled                    
            wnbpts=size(clustens{dispUnitSel{1}(1)},3);
            t=(1:wnbpts)';                    
            ts=linspace(1,wnbpts,wnbpts*obj.dwsamplefactor)';
            xdata=linspace(0,wnbpts*obj.dwsamplefactor/(obj.Params.fileinfo.samplingrate),wnbpts*obj.dwsamplefactor)'*1000;%size(clustens,4));
        else
            xdata=linspace(0,size(clustens{dispUnitSel{1}(1)},3)*obj.dwsamplefactor/(obj.Params.fileinfo.samplingrate),size(clustens{dispUnitSel{1}(1)},3))'*1000;%size(clustens,4));
        end
        for u=1:nbunit
            hold(hwaveforms(1,u),'off');
            if Frealign && nbunit==2
                [crossval ttu]=obj.MultiChSpkCrossCorr(UnitSel);
                xdataU=xdata(ttu{1,2}(1,:));
            else
                xdataU=xdata;
            end 
            for c=1:size(dispUnitSel{u},2)
                if dispUnitSel{u}(c)>0 && dispUnitSel{u}(c)<=size(clustens,1)
                    optimCh=round(obj.MclusterChID(dispUnitSel{u}(c)));
                    idx=(find(max(squeeze(clustens{dispUnitSel{u}(c)}(optimCh,:,:))')~=0))';
                    if numel(idx)==1
                        idx=[idx idx];
                    end
                    if ~Fwhitened
                        clustens{dispUnitSel{u}(c)}=single(clustens{dispUnitSel{u}(c)})/obj.Params.detect.SpkADgain;
                    else
                        %whitenclust=clustens{dispUnitSel{u}(c)}(:,idx,:);
                        whitenclust=obj.NoiseWhitening(clustens{dispUnitSel{u}(c)}(:,idx,:));
                    end
%                             whitenclusttemp=zeros(253,size(whitenclust,2),189);
%                             for k=1:size(whitenclust,2)
%                                 whitenclusttemp(:,k,:)=interp2(squeeze(whitenclust(:,k,:)),2);
%                             end
%                             whitenclust=whitenclusttemp;
                    maxclust=0;                            
                    for ch=1:nbCh                        
                        ydata=[];                                                        
                        %idy=(find(max(squeeze(clustens(ChSel(ch),UnitSel(u),:,:)))~=0))';
                        if ~isempty(idx)
                            if Frealign && nbunit==2
                                if Fwhitened
                                    ydata=[ydata;squeeze(whitenclust(ChSel(ch),:,:))];
                                else
                                    ydata=squeeze(clustens{dispUnitSel{u}(c)}(ChSel(ch),idx,:));
                                end
                                ydata=ydata(:,ttu{1,2}(u,:));
                                xdataU=xdata(ttu{1,2}(1,:));
                            elseif Foversampled
                                if Fwhitened
                                    x = squeeze(whitenclust(ChSel(ch),:,:))';
                                else
                                    x = squeeze(clustens{dispUnitSel{u}(c)}(ChSel(ch),idx,:))';
                                end
                                ydata= (sinc(ts(:,ones(size(t))) - t(:,ones(size(ts)))')*x)';                                        
                            else
                                if Fwhitened
                                    ydata=[ydata;squeeze(whitenclust(ChSel(ch),:,:))];
                                else
                                    ydata=[ydata;squeeze(clustens{dispUnitSel{u}(c)}(ChSel(ch),idx,:))];
                                end                                        
                            end
                        end
                        if ~isempty(ydata)
                            stdydata=std(ydata);
                            meanydata=mean(ydata);
                            if max(max(abs(meanydata)))>maxclust
                                maxclust=max(max(abs(meanydata)));
                            end
%                                     for k=1:size(ydata,1)
%                                         Pxx = fft(ydata(k,:));
%                                         ydata(k,:)=abs(Pxx);
%                                     end
                            %xdataU=w/(2*pi)*4000;  

%                                     vf0=FiltFiltM(butHB,butHA,ydata',1);
%                                     filydata=(FiltFiltM(butLB,butLA,vf0))';
%                                     ydata=filydata;

                            maxclust=max(maxclust,max(abs(meanydata)));
                            if Fdispstd
                                fill([xdataU' fliplr(xdataU')],[meanydata-stdydata-maxY*(ch-1) fliplr(meanydata+stdydata-maxY*(ch-1))],[0.8 0.8 0.8],'Parent',hwaveforms(1,u),'LineStyle','none'); 
                            else
%                                         xdataU=linspace(0,6,size(ydata,2));
                                plot(hwaveforms(1,u),xdataU,ydata'-maxY*(ch-1),'Color',[0.8 0.8 0.8]);
                            end
                            set(hwaveforms(1,u),'NextPlot','add');
                            plot(hwaveforms(1,u),xdataU,meanydata-maxY*(ch-1),'Color',[1 0 0],'LineWidth',1.5,'LineStyle','--');                                    
                        end
                        hold(hwaveforms(1,u),'on');
                    end
                    %disp(maxclust)
                    %fprintf([num2str(maxclust) ' ']);
                end                        
            end 
            set(hwaveforms(1,u),'YLim',[-maxY*nbCh maxY]);
            set(hwaveforms(1,u),'XLim',[min(xdataU) max(xdataU)]);
            if u==1 
                set(hwaveforms(1,u),'YTick',-maxY*(nbCh-1):maxY:maxY);
                set(hwaveforms(1,u),'YTickLabel',num2str([fliplr(ChSel) maxY]'));
            else
                set(hwaveforms(1,u),'YTick',[]);
            end
            try
                set(get(hwaveforms(1,u),'Title'),'String',['clust#' num2str(UnitSel(u)) ' SU#' num2str(obj.getUnitIDMultiCh(UnitSel(u)))...
                                                            ' Ch#' num2str(obj.MclusterChID(UnitSel(u)),2+(obj.MclusterChID(UnitSel(u))>9))],'FontSize',8);
                set(get(hwaveforms(1,u),'XLabel'),'String','ms');
            catch
                warning('no more clusters')
            end
        end
%                 f=figure;set(hwaveforms(1,u),'Parent',f)
%                 savefig(f,['/storage/laur/Data_2/FournierJulien/SpikeSorting ground truth/hc-1/new-hc1/figures/final/K_27_1/clusters/cluster' num2str(UnitSel(u)) '.fig']);
    end            
end