function obj = Convert2NewFileFormat(obj,filepath)
    if nargin<2
        if ischar(obj.filepath)
            filepath=uigetdir(obj.filepath);
        else
            filepath=uigetdir(['D:' filesep]);
        end
    end
    nameposit=strfind(filepath,filesep);
    if ~isempty(nameposit)
        k=length(nameposit);
        filname=filepath(nameposit(k)+1:size(filepath,2));
    end

    obj.filepath=filepath;
    obj.filename=filname; 
    obj.Params.fileinfo=[];
        
    obj.LoadSpikeData();
    obj.LoadSortedSpikes();
    
    infoNew.fileinfo.acqsystem=obj.acqsystem;
    infoNew.fileinfo.samplingrate=obj.samplingrate;
    infoNew.fileinfo.ADgain=obj.ADgain;
    infoNew.fileinfo.H0gain=obj.H0gain;
    infoNew.fileinfo.InputRange=obj.InputRange;
    infoNew.fileinfo.Channelcount=obj.Channelcount;
    infoNew.fileinfo.Channelnum=obj.Channelnum;
    infoNew.fileinfo.electrode=obj.Params.detect.electrode;
    infoNew.fileinfo.ChannelPosX=obj.ChannelPosX;
    infoNew.fileinfo.ChannelPosY=obj.ChannelPosY;
    load([obj.filepath filesep obj.filename],'fileinfo');
    
    fileinfo.InputRange = infoNew.fileinfo.InputRange;
    fileinfo.electrode = infoNew.fileinfo.electrode;
    fileinfo.ChannelPosX = infoNew.fileinfo.ChannelPosX;
    fileinfo.ChannelPosY = infoNew.fileinfo.ChannelPosY;
    save([obj.filepath filesep obj.filename '.mat'],'fileinfo','-append');
    
    
    infoNew.detect.length2load=obj.Params.detect.length2load;
    infoNew.detect.FcHigh=obj.Params.detect.FcHigh;
    infoNew.detect.FcLow=obj.Params.detect.FcLow;
    infoNew.detect.butterOrder=obj.butterOrder;
    infoNew.detect.Fsubmedian=obj.Params.detect.Fsubmedian;
    infoNew.detect.Fspline=obj.Params.detect.Fspline;
    infoNew.detect.zthresh=obj.Params.detect.zthresh;
    infoNew.detect.wtime1=obj.Params.detect.wtime1;
    infoNew.detect.wtime2=obj.Params.detect.wtime2;
    infoNew.detect.SpkADgain=obj.Params.detect.SpkADgain;
    infoNew.detect.Tcensoredfactor=obj.Params.detect.Tcensoredfactor;
    infoNew.detect.FGroundtruth=obj.Params.detect.FGroundtruth;

    
    infoNew.MCluster.TrainingMaxnbSpk=obj.Params.MCluster.TrainingMaxnbSpk;
    infoNew.MCluster.nbPC=obj.Params.MCluster.nbPC;
    infoNew.MCluster.UpdateTime=obj.Params.MCluster.UpdateTime;
    infoNew.MCluster.ProjRange=obj.Params.MCluster.ProjRange;
    infoNew.MCluster.NbIteration = 200;
    infoNew.MCluster.nbmaxunits = floor(50000^0.5);
    infoNew.MCluster.minUsize = obj.Params.MCluster.minUsize;
    infoNew.MCluster.bestminUsize=obj.Params.MCluster.minUsize;
    infoNew.MCluster.Pmisthreshold=obj.Params.MCluster.Pmisthreshold;
    infoNew.MCluster.Fgreedy=obj.Params.MCluster.Fgreedy;
    infoNew.MCluster.cdfChi2thresh=obj.Params.MCluster.cdfChi2thresh;
    
    obj.Params = infoNew;
    
    spkwave=obj.Spkwave;
    nbCh=size(spkwave,1);
    nbunitstotal = numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0)));
    clust1=cell(nbCh,nbunitstotal);
    for u=1:nbunitstotal
        spkevtidx=find(obj.Spkevent(6,:)==u & obj.Spkevent(10,:)==0);
        for ch=1:nbCh
            if numel(spkevtidx)>1
                clust1{ch,u}=mean(squeeze(spkwave(ch,spkevtidx,:)));
            else
                clust1{ch,u}=(spkwave(ch,spkevtidx,:));
            end
        end
    end
    obj.SpkwaveclustAve = cell(1,nbunitstotal);
    for u=1:nbunitstotal
        obj.SpkwaveclustAve{ordposidx(u)} = cell2mat(clust1(:,u));
    end
    
    
    
    fname=[obj.datafolder filesep obj.filename '_Spkdata.mat'];
    info=obj.Params;
    spkevent=obj.Spkevent;
    spkwave=obj.Spkwave;
    spkmask=obj.Spkmask;

    whiteningmatallspatial=obj.WhiteningMatAllSpatial;
    rawdata=int16(obj.RawData(:,1:floor(0.5*end))*obj.Params.detect.SpkADgain);
    save(fname,'spkevent','spkwave','spkmask','whiteningmatallspatial','rawdata','info','-v7.3');
    
    fname=[obj.datafolder filesep obj.filename '_Spkevent.mat'];
    info=obj.Params;
    spkevent=obj.Spkevent;
    spkeventUM=obj.SpkeventUM;
    idxglobalAll=obj.IDXglobalAll;
    errorglobalAll=obj.ErrorglobalAll;
    quality=obj.Quality;
    spkclustIDmulti=obj.SpkclustIDMulti;
    spkfeature=obj.SpkFeatures;
    
    munitchID=obj.MunitChID;
    mclusterchID=obj.MclusterChID;
    nbunitmultich=obj.nbUnitMultiCh;
    spkwaveclustave=obj.SpkwaveclustAve;
        
    save(fname,'spkevent','spkeventUM','idxglobalAll','errorglobalAll','quality','spkclustIDmulti','spkwaveclustave','munitchID','spkfeature','mclusterchID','nbunitmultich','info','fcomment','-v7.3');
end