function obj=LoadDATFile(obj,num)
%Loads .dat file data between timestamps obj.Epstartime(num) and
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

    fname=[obj.filepath filesep obj.filename '.dat'];    
    
    if strcmp(obj.filepath,'')~=1
        datfile = fopen(fname, 'r');
        fseek(datfile,(obj.timestamp1-1)*obj.Params.fileinfo.Channelcount*2,'bof');
        obj.RawData = double(fread(datfile,[obj.Params.fileinfo.Channelcount obj.timestamp2-obj.timestamp1],'int16'));
        obj.RawData = obj.RawData(obj.Params.fileinfo.Channelnum,:);
        if obj.timestamp2 >= size(obj.RawData,2)
            obj.timestamp2 = size(obj.RawData,2);
        end
        fclose(datfile);
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