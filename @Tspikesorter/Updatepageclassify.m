function obj=Updatepageclassify(obj,UnitSel,Fwhitened)
    parentpanel=obj.pageclassify;
    if Fwhitened
        maxY=15;
    else
        maxY=80;
    end
    nbunit=length(UnitSel);

    Clustersubset=cell(1,length(UnitSel));
    unitChID=zeros(1,length(UnitSel));
    for u=1:length(UnitSel)
        Clustersubset{u}=obj.getClusterIDMultiCh(UnitSel(u));
    end
    if ismember(300,UnitSel)
        Clustersubset{UnitSel==300}=-1;  
    end

    ChSel=[];
    for u=1:nbunit
        for c=1:size(Clustersubset{u},2)
            ChSel=[ChSel (round(obj.MclusterChID(Clustersubset{u}(c)))-4):(round(obj.MclusterChID(Clustersubset{u}(c)))+4)];
        end
    end
    ChSel=unique(ChSel(ChSel>1 & ChSel<size(obj.Spkwave,1)));
    nbCh=size(ChSel,2);

    spkevtU=obj.SpkeventUM;            

    hwaveforms=zeros(nbunit,4);            

    if ~isempty(get(parentpanel,'Children'))
        delete(get(parentpanel,'Children'));
    end
    for u=1:nbunit
        hwaveforms(u,1)=subplot(3,2*nbunit,[1 2*nbunit+1 2*nbunit*2+1]+(u-1)*2,'Parent',parentpanel);
        hwaveforms(u,2)=subplot(3,2*nbunit,u*2,'Parent',parentpanel);
        hwaveforms(u,3)=subplot(3,2*nbunit,u*2+2*nbunit,'Parent',parentpanel);
        hwaveforms(u,4)=subplot(3,2*nbunit,u*2+2*nbunit*2,'Parent',parentpanel);                    
    end

    spktime=cell(1,nbunit);
    spkamp=cell(1,nbunit);
    varexpl=zeros(1,nbunit);
    spktimend=(spkevtU(2,end)-spkevtU(2,1))*10^-3;

    for u=1:nbunit
        spktimeU=find(ismember(spkevtU(6,:),Clustersubset{u}));
        if ~isempty(spktimeU)     
            fileID=unique(spkevtU(1,:));
            spktime{u}=(spkevtU(2,spkevtU(1,:)==fileID(1) & ismember(spkevtU(6,:),Clustersubset{u}))-spkevtU(2,1))*10^-3; 
            spkamp{u}=spkevtU(5,spktimeU);
            offset=0;
            if size(fileID,2)>1
                for k=2:size(fileID,2)
                    if ~isempty(spktime{u})
                        offset=offset+spkevtU(2,find(spkevtU(1,:)==fileID(k-1),1,'last'))*10^-3;
                    end 
                    spktime{u}=[spktime{u} offset+(spkevtU(2,spkevtU(1,:)==fileID(k) & ismember(spkevtU(6,:),Clustersubset{u}))-spkevtU(2,1))*10^-3];                            
                end
            end                        
            if ~isempty(spktime{u})
                xisi=(0:0.5:10);
                if length(spktime{u})>1
                    spkisi=histc(diff(spktime{u}),xisi);
                else
                    spkisi=histc(10000,xisi);
                end
                xpsth=(0:1000:(spktimend))*10^-3;
                spkpsth=histc(spktime{u}*10^-3,xpsth);
                xcrosscorrU=-50:0.5:50;
                crosscorrU=zeros(1,length(xcrosscorrU));
                spkvec1=round(spktime{1})-0.5*((round(spktime{1})-round(spktime{1}-0.25))>0);
                spkvecu=round(spktime{u})-0.5*((round(spktime{u})-round(spktime{u}-0.25))>0);
                for tt=1:length(xcrosscorrU)
                    crosscorrU(tt)=(xcrosscorrU(tt)~=0)*sum(ismember(spkvecu+xcrosscorrU(tt),spkvec1));
                end
                varexpl(u)=(1-mean(spkevtU(8,ismember(spkevtU(6,:),Clustersubset{u})))/mean(spkevtU(9,ismember(spkevtU(6,:),Clustersubset{u}))))*100;
            end
            bar(hwaveforms(u,2),xisi,spkisi);
            scatter(hwaveforms(u,3),spktime{u},spkamp{u});
            bar(hwaveforms(u,4),xcrosscorrU,crosscorrU);

%                     f=figure;
%                     bar(xcrosscorrU,crosscorrU);
%                     set(gca,'XLim',[-50 50]);
%                     savefig(f,['/storage/laur/Data_2/FournierJulien/SpikeSorting ground truth/hc-1/figure/ctrl_consensus_cluster' num2str(Clustersubset{u}) '_ISI.fig']);

            set(hwaveforms(u,2),'XLim',[0 10]);
            set(hwaveforms(u,4),'XLim',[-50 50]);
            set(hwaveforms(u,3),'XLim',[-200 max(spktime{u})+200]);

            set(get(hwaveforms(u,2),'Title'),'String',['ISI' ' "clust# ' num2str(Clustersubset{u}(1)) '-' num2str(Clustersubset{u}(end)) '"'],'FontSize',8);
            set(get(hwaveforms(u,3),'Title'),'String',['PSTH' ' unit#' num2str(UnitSel(u)) ' clust#' num2str(Clustersubset{u})...
                                                                  ' optCh#' num2str(obj.MclusterChID(max(Clustersubset{u},1)))],'FontSize',8);
            set(get(hwaveforms(u,4),'Title'),'String',['CrossCorr' ' "' num2str(Clustersubset{u}(1)) '-' num2str(Clustersubset{u}(end))  '"vs"' num2str(UnitSel(find(unitChID==1,1,'first'))) '-' num2str(UnitSel(find(unitChID==1,1,'last'))) '"'],'FontSize',8);
            set(get(hwaveforms(u,2),'XLabel'),'String','ms');
            set(get(hwaveforms(u,3),'XLabel'),'String','sec');
            set(get(hwaveforms(u,4),'XLabel'),'String','ms');

            xdata=linspace(0,size(obj.Spkwave,3)/(obj.Params.fileinfo.samplingrate/obj.dwsamplefactor),size(obj.Spkwave,3))'*1000;%size(clustens,4));

            hold(hwaveforms(u,1),'off');
            idxspkoverlaps=find(diff(spktime{u})<2);
            if isempty(idxspkoverlaps)
                idxspkoverlaps=1:numel(spktime{u});
            end
            idxspku=find(ismember(obj.Spkevent(6,:),Clustersubset{u}));
            idxspku=idxspku(1:min(end,300));
            if ~Fwhitened
                clustens=single(obj.Spkwave(:,idxspku,:))/obj.Params.detect.SpkADgain;
            else
                clustens=obj.NoiseWhitening(obj.Spkwave(:,idxspku,:));%([1 idxspkoverlaps]);
                %clustens=clustens(:,53:65,:);
                clustens(:,1,:)=[];
            end
            maxclust=0;
            for ch=1:nbCh                       
                ydata=squeeze(clustens(ChSel(ch),:,:));                        
                if ~isempty(ydata)
                    stdydata=std(ydata);
                    meanydata=mean(ydata);
                    if max(max(abs(meanydata)))>maxclust
                        maxclust=max(max(abs(meanydata)));
                    end
                    plot(hwaveforms(u,1),xdata,ydata-maxY*(ch-1),'Color',[0.8 0.8 0.8]);
                    %fill([xdata' fliplr(xdata')],[meanydata-stdydata-maxY*(ch-1) fliplr(meanydata+stdydata-maxY*(ch-1))],[0.8 0.8 0.8],'Parent',hwaveforms(u,1),'LineStyle','none');
                    set(hwaveforms(u,1),'NextPlot','add');
                    plot(hwaveforms(u,1),xdata,meanydata-maxY*(ch-1),'Color',[1 0 0],'LineWidth',1.5,'LineStyle','--');
                end
                hold(hwaveforms(u,1),'on');
            end

            set(hwaveforms(u,1),'YLim',[-maxY*nbCh maxY]);
            set(hwaveforms(u,1),'XLim',[-2 8]);%min(xdata) max(xdata)]);
            if u==1
                set(hwaveforms(u,1),'YTick',-maxY*(nbCh-1):maxY:maxY);
%                         set(hwaveforms(u,1),'YTickLabel',num2str([wrev(ChSel) maxY]'));
            else
                set(hwaveforms(u,1),'YTick',[]);
            end
            set(get(hwaveforms(u,1),'XLabel'),'String','ms');                    
        end                
    end                                                       
end