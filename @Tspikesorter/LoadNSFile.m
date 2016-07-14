function obj=LoadNSFile(obj,num)
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

    fname=[obj.filepath filesep obj.filename '.ns5'];    
    
    global pepNEV;
    if strcmp(obj.filepath,'')~=1   
        if isempty(pepNEV)
            [tmp,SamplingRateInKHZ,nchan] = nsopen(fname);
        end
        if obj.timestamp2 <= size(pepNEV.ns.Data.data,2)
            obj.RawData = double(pepNEV.ns.Data.data(obj.Params.fileinfo.Channelnum,obj.timestamp1:obj.timestamp2));
        else
            obj.RawData = double(pepNEV.ns.Data.data(obj.Params.fileinfo.Channelnum,obj.timestamp1:end));
            obj.timestamp2 = size(pepNEV.ns.Data.data,2);
        end
    end
        
    %timestamps are in sample units for MCD files. We convert it in
    %microseconds after loading the data as later we'll use microseconds
    %throughout
    obj.timestamp1 = obj.timestamp1 / (obj.Params.fileinfo.samplingrate)*10^6; %obj.Epstartime(num)/(obj.Params.fileinfo.samplingrate)*10^6;   
    obj.timestamp2 = obj.timestamp2 / (obj.Params.fileinfo.samplingrate)*10^6; %
%     if num+1<=numel(obj.Epstartime)
%         obj.timestamp2=obj.Epstartime(num+1)/(obj.Params.fileinfo.samplingrate)*10^6;
%     else
%         if num~=1
%             obj.timestamp2=(obj.Epstartime(num)+(obj.Epstartime(num)-obj.Epstartime(num-1)))/(obj.Params.fileinfo.samplingrate)*10^6;
%         else
%             obj.timestamp2=size(obj.RawData,2)/(obj.Params.fileinfo.samplingrate)*10^6;
%         end
%     end    
end