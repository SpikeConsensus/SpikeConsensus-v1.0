function obj=FilterRawData(obj)
% band pass filter the raw data contained in obj.RawData between 
% Params.detect.Fclow and Params.detect.Fchigh
% If obj.Params.detect.Fsubmedian is true, we also do a median filtering on
% each probe.
    Saturation=obj.Params.fileinfo.InputRange;
    idxSat=cell(1,obj.Params.fileinfo.Channelcount);
              
    vf0=filtfilt(obj.butterHB,obj.butterHA,obj.RawData');
    obj.RawData=(filtfilt(obj.butterLB,obj.butterLA,vf0))';
    
    
    if obj.Params.detect.Fsubmedian
        %copy of obj.RawData for parfor loops   
        rawdata=obj.RawData;  
        %we check if any channel saturated. If so, we'll
        %ignore this segment of data when doing the median filter 
        for k=1:obj.Params.fileinfo.Channelcount
            idxSat{k}=find(rawdata(k,:)>=0.999*Saturation & rawdata(k,:)<=1.001*Saturation);
            if ~isempty(idxSat{k})
                disp(['saturation detected on Channel#' num2str(k)]);
            end
        end
        rawdataNorm=zeros(size(rawdata));
        for el=1:size(obj.elconfig,2)
            parfor ch=obj.elconfig(1,el):obj.elconfig(2,el) 
                rawdataNorm(ch,:)=rawdata(ch,:)/(median(abs(rawdata(ch,:))/0.6745));
            end
            VmedianCh=median(rawdataNorm(obj.elconfig(1,el):obj.elconfig(2,el),:));
            VmedianChNorm=sum(VmedianCh.^2);

            parfor ch=obj.elconfig(1,el):obj.elconfig(2,el) 
                rawdata(ch,:)=rawdata(ch,:)-sum(rawdata(ch,:).*VmedianCh)*VmedianCh/(VmedianChNorm);            
            end
            for ch=obj.elconfig(1,el):obj.elconfig(2,el)
                if ~isempty(idxSat{ch})
                    rawdata(obj.elconfig(1,el):obj.elconfig(2,el),idxSat{ch})=0;
                end
            end
        end
        rawdata=(FiltFiltM(obj.butterLB,obj.butterLA,rawdata'))';
        obj.RawData=rawdata;
    end        
end