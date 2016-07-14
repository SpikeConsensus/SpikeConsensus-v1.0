function obj=LoadNeuralynxCSC(obj,num)
%Loads Neuralynx CSC file data between timestamps obj.Epstartime(num) and
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
    fname=cell(1,obj.Params.fileinfo.Channelcount);
    loadParams=cell(1,obj.Params.fileinfo.Channelcount);
    for k=1:obj.Params.fileinfo.Channelcount
        ch_num=obj.Params.fileinfo.Channelnum(k);
        ch_str=int2str(ch_num);
        fname{k}=[obj.datafolder filesep 'CSC' ch_str '.ncs'];
        loadParams{k}=struct;
        loadParams{k}.FieldSelectionFlags= [1 0 1 0 1]; % [TimeStamps, ChannelNumbers, SampleFreq, NumberOfValidSamples, Samples];
        loadParams{k}.HeaderExtractionFlag = 0;
        loadParams{k}.ExtractionMode = 4;
        loadParams{k}.ExtractionModeVector = [obj.timestamp1 obj.timestamp2];
        loadParams{k}.samplingrate=obj.Params.fileinfo.samplingrate;
    end

    Header= NlxCSC2mat(fname{1},[0 0 0 0 0],1,1,[]);
    Linputrange=1;
    while strcmp(Header{Linputrange}(1:min(end,11)),'-InputRange')~=1
        Linputrange=Linputrange+1;
    end
    obj.Params.fileinfo.InputRange=str2double(Header{Linputrange}(13:end));

    if ~(obj.Params.fileinfo.InputRange>=1000 && obj.Params.fileinfo.InputRange<=3000)
        error('bad input range detected')
    end

    for k=1:obj.Params.fileinfo.Channelcount
        if strcmp(obj.filepath,'')~=1
            if ~exist(fname{k})
                warning(['file does not exist: ' fname{k}]);
            else
                [TimeStamps,SamplingFq,Samples]=NlxCSC2mat(fname{k},loadParams{k}.FieldSelectionFlags,loadParams{k}.HeaderExtractionFlag,loadParams{k}.ExtractionMode, loadParams{k}.ExtractionModeVector);
                startindex=floor((loadParams{k}.ExtractionModeVector(1)-TimeStamps(1))*10^(-6)*loadParams{k}.samplingrate);
                endindex=512-floor((loadParams{k}.ExtractionModeVector(2)-TimeStamps(end))*10^(-6)*loadParams{k}.samplingrate);
                dataRaw{k} = reshape(Samples',1,size(Samples,1)*size(Samples,2));
                if endindex>0
                    dataRaw{k}((end-endindex):end)=[];
                elseif endindex<0
                    %dataRaw{k}=[dataRaw{k} zeros(size(dataRaw{k},1),(-endindex))];
                end
                if startindex>0
                    dataRaw{k}(1:startindex-1)=[];
                elseif startindex<0
                    dataRaw{k}=[zeros(size(dataRaw{k},1),-startindex) dataRaw{k}];
                end
            end
            if obj.Params.detect.Fspline
                nbinterpfactor=4;
                nbptsinit=size(dataRaw{k},2);
                dataRaw{k}=spline(1:nbptsinit,dataRaw{k},linspace(1,nbptsinit,nbinterpfactor*nbptsinit));                
            end
        end
    end    
    obj.RawData=double(cell2mat(dataRaw)*scaleFactor*1e6/gain);
end

