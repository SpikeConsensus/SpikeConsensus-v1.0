function obj=Updatepagedetect(obj,ChSel)
    parentpanel=obj.pagedetect;
    nbCh=size(ChSel,2);
    hrawdata=zeros(1,nbCh);
    hrawwaveforms=zeros(1,nbCh);
    hwaveAmplitudes=zeros(1,nbCh);

    if ~isempty(get(parentpanel,'Children'))
        delete(get(parentpanel,'Children'));
    end
    for ch=1:nbCh
        hrawdata(ch)=subplot(2,2*nbCh,[(ch-1)*2+1 (ch-1)*2+2],'Parent',parentpanel);
        hrawwaveforms(ch)=subplot(2,2*nbCh,nbCh*2+(ch-1)*2+1,'Parent',parentpanel);
        hwaveAmplitudes(ch)=subplot(2,2*nbCh,nbCh*2+(ch-1)*2+2,'Parent',parentpanel);                
    end
    if ~isempty(obj.RawData)
        xdata=linspace(0,size(obj.RawData,2)/obj.Params.fileinfo.samplingrate,size(obj.RawData,2));
        for ch=1:nbCh
            if ChSel(ch)<=size(obj.RawData,1)
                plot(hrawdata(ch),xdata,single(obj.RawData(ChSel(ch),:))/obj.Params.detect.SpkADgain,'Color',[0 0 0]);
                set(hrawdata(ch),'NextPlot','add');
%                 plot(hrawdata(ch),xdata,obj.detectthresh(ChSel(ch))*ones(size(xdata)),'Color','r');
%                 plot(hrawdata(ch),xdata,-obj.detectthresh(ChSel(ch))*ones(size(xdata)),'Color','r');
                set(get(hrawdata(ch),'Title'),'String',['Raw data ch#' num2str(ChSel(ch))]); 
            end
        end
    end            
    if ~isempty(obj.Spkwave)
        xdata=linspace(0,size(obj.Spkwave,3)/(obj.Params.fileinfo.samplingrate/obj.dwsamplefactor),size(obj.Spkwave,3));
        for ch=1:nbCh
            if ChSel(ch)<=size(obj.Spkwave,1)
                idx=find(obj.Spkevent(3,:)==ChSel(ch),1000,'first');
                idx(idx>size(obj.Spkwave,2))=[];
                maxval=max(max(single(squeeze(obj.Spkwave(ChSel(ch),idx,:)))/obj.Params.detect.SpkADgain));
                %plot(hrawwaveforms(ch),xdata,squeeze(obj.Spkwave(ChSel(ch),idx(obj.Spkwave(ChSel(ch),idx,17)==maxval),:))');
                if ~isempty(idx)
                    plot(hrawwaveforms(ch),xdata,single(squeeze(obj.Spkwave(ChSel(ch),idx,:))')/obj.Params.detect.SpkADgain);
                    set(hrawwaveforms(ch),'NextPlot','add');
                    plot(hrawwaveforms(ch),xdata,std(single(squeeze(obj.Spkwave(ChSel(ch),idx,:))))/obj.Params.detect.SpkADgain,'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
                    set(get(hrawwaveforms(ch),'Title'),'String',['waveforms ch#' num2str(ChSel(ch))]);
                else
                    plot(hrawwaveforms(ch),xdata,zeros(size(xdata))/obj.Params.detect.SpkADgain);
                    set(hrawwaveforms(ch),'NextPlot','add');
                    plot(hrawwaveforms(ch),xdata,zeros(size(xdata))/obj.Params.detect.SpkADgain,'Color',[0 0 0],'LineWidth',2,'LineStyle','--');
                    set(get(hrawwaveforms(ch),'Title'),'String',['waveforms ch#' num2str(ChSel(ch))]);
                end
            end
        end                
    end
    if ~isempty(obj.Spkevent)
        for ch=1:nbCh
            if ChSel(ch)<=size(obj.Spkwave,1)
                xdata=(obj.Spkevent(2,obj.Spkevent(3,:)==ChSel(ch))-obj.Spkevent(2,1))*10^-6;
                ydata=obj.Spkevent(5,obj.Spkevent(3,:)==ChSel(ch));
                plot(hwaveAmplitudes(ch),xdata,ydata,'LineStyle','none','Marker','+');
                set(get(hwaveAmplitudes(ch),'Title'),'String',['SNR (z-score) ch#' num2str(ChSel(ch))]);
            end
        end                
    end                     
end