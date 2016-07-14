function obj=LoadMcdFile(obj,num)
%Loads MCD file data between timestamps obj.Epstartime(num) and
%obj.Epstartime(num+1) into obj.RawData

    scaleFactor= obj.Params.fileinfo.ADgain;
    gain = obj.Params.fileinfo.H0gain;
    
    obj.timestamp1=obj.Epstartime(num);
    
    if num+1<=numel(obj.Epstartime)
        obj.timestamp2=obj.Epstartime(num+1);
    else
        obj.timestamp2=obj.Epstartime(num)+(obj.Epstartime(num)-obj.Epstartime(num-1));
    end

    dataRaw=cell(obj.Params.fileinfo.Channelcount,1);

    fname=[obj.filepath filesep obj.filename '.mcd'];    
    
    [~, hfile] = ns_OpenFile(fname);
    [~,info] = ns_GetFileInfo(hfile); 
    nbSamples = info.TimeSpan/info.TimeStampResolution;
    
    if strcmp(obj.filepath,'')~=1
        for k=1:obj.Params.fileinfo.Channelcount   
            [~,~,data]=ns_GetAnalogData(hfile,obj.Params.fileinfo.Channelnum(k),floor(obj.timestamp1),min(nbSamples-floor(obj.timestamp1)+1,floor(obj.timestamp2-obj.timestamp1)));
            dataRaw{k}=data';   
            if obj.Params.detect.Fspline
                nbinterpfactor=4;
                nbptsinit=size(dataRaw{k},2);
                dataRaw{k}=spline(1:nbptsinit,dataRaw{k},linspace(1,nbptsinit,nbinterpfactor*nbptsinit));                
            end
        end
    end    
    obj.RawData=double(cell2mat(dataRaw)*scaleFactor*10^6/gain);
    
    if obj.Params.detect.FGroundtruth
        load([fname(1:(find(fname=='.',1,'last')-1)) '_intraAP.mat'],'intraAP');
        obj.Vmintra=double(intraAP(floor(obj.timestamp1):min(numel(intraAP),floor(obj.timestamp2-1))));
    end
    
    %timestamps are in sample units for MCD files. We convert it in
    %microseconds after loading the data as later we'll use microseconds
    %throughout
    obj.timestamp1=obj.Epstartime(num)/(obj.Params.fileinfo.samplingrate)*10^6;    
    if num+1<=numel(obj.Epstartime)
        obj.timestamp2=obj.Epstartime(num+1)/(obj.Params.fileinfo.samplingrate)*10^6;
    else
        if num~=1
            obj.timestamp2=(obj.Epstartime(num)+(obj.Epstartime(num)-obj.Epstartime(num-1)))/(obj.Params.fileinfo.samplingrate)*10^6;
        else
            obj.timestamp2=size(obj.RawData,2)/(obj.Params.fileinfo.samplingrate)*10^6;
        end
    end    
end