function obj=LoadMatFile(obj,num)
%Loads Mat file data between timestamps obj.Epstartime(num) and
%obj.Epstartime(num+1) into obj.RawData.
    stamp1=max(1,round(obj.Epstartime(2,num)));    
    if num+1<=size(obj.Epstartime,2)
        stamp2=obj.Epstartime(2,num+1);
    else
        if num~=1
            stamp2=obj.Epstartime(2,num)+(obj.Epstartime(2,num)-obj.Epstartime(2,num-1));%obj.datafilelength(num)*10^6+obj.sidetime*10^6;
        else
            stamp2=+inf;
        end
    end
    stamp2=round(stamp2);
    load([ obj.filepath filesep obj.filename '_data.mat'],'data');
    stamp2=min(size(data,2),stamp2);
    obj.RawData=double(data(1:obj.Params.fileinfo.Channelcount,stamp1:stamp2))/10;     
    
    if obj.Params.detect.FGroundtruth
        ChannelIntra=obj.Params.fileinfo.Channelcount+1:size(data,1);
        if abs(min(min(data(ChannelIntra,:))))>200
            Vmgain=10;
        else
            Vmgain=1;
        end
        obj.Vmintra=double(data(ChannelIntra,stamp1:min(size(data,2),stamp2)))/Vmgain;
    end
        
    if obj.Params.detect.Fspline
        nbinterpfactor=4;
        rawdata=obj.RawData;
        nbptsinit=size(rawdata,2);
        obj.RawData=zeros(size(rawdata,1),nbinterpfactor*nbptsinit);
        for ch=1:size(obj.RawData,1)
            obj.RawData(ch,:)=spline(1:nbptsinit,rawdata(ch,:),linspace(1,nbptsinit,nbinterpfactor*nbptsinit));
        end
        vmintra=obj.Vmintra;
        nbptsinit=size(vmintra,2);
        obj.Vmintra=zeros(size(vmintra,1),nbinterpfactor*nbptsinit);
        for u=1:size(obj.Vmintra,1)
            obj.Vmintra(u,:)=spline(1:nbptsinit,vmintra(u,:),linspace(1,nbptsinit,nbinterpfactor*nbptsinit)); 
        end
    else
        nbinterpfactor=1;
    end
    
    obj.timestamp1=obj.Epstartime(2,num)/(obj.Params.fileinfo.samplingrate/nbinterpfactor)*10^6;    
    if num+1<=size(obj.Epstartime,2)
        obj.timestamp2=obj.Epstartime(2,num+1)/(obj.Params.fileinfo.samplingrate/nbinterpfactor)*10^6;
    else
        if num~=1
            obj.timestamp2=(obj.Epstartime(2,num)+(obj.Epstartime(2,num)-obj.Epstartime(2,num-1)))/(obj.Params.fileinfo.samplingrate/nbinterpfactor)*10^6;%obj.datafilelength(num)*10^6+obj.sidetime*10^6;
        else
            obj.timestamp2=size(obj.RawData,2)/(obj.Params.fileinfo.samplingrate/nbinterpfactor)*10^6;
        end
    end
    obj.Params.fileinfo.InputRange = +inf;
end

