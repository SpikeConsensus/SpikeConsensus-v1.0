function obj=SpkClusterMultiCh(obj,FKitespace)
% performs N iterations of Kmeans clustering and saves in IDXglobalAll 
% and ErrorglobalAll the cluster labels and the associated chi2 values 
% of each spike, for every iteration.

    if nargin < 2
        FKitespace = false;
    end
    
    wnbptsinit=size(obj.Spkwave,3);
    
    if obj.Params.detect.Tcensoredfactor<=0.25
        obj.Params.MCluster.Timewindowfactor = 3*obj.Params.detect.Tcensoredfactor;
    else
        obj.Params.MCluster.Timewindowfactor = 2*obj.Params.detect.Tcensoredfactor;
    end
    wnbpts=floor(wnbptsinit*obj.Params.MCluster.Timewindowfactor/2)*2+1;
    wnbpts=min(wnbpts,wnbptsinit);  
    timeresol=wnbpts/2*obj.dwsamplefactor/obj.Params.fileinfo.samplingrate*10^6;
        
    obj.SpkFeatures=zeros(obj.Params.MCluster.nbPCmax*size(obj.Spkwave,1),size(obj.Spkwave,1),wnbpts);
    totalPC=0;            
    obj.Spkevent(6,:)=0*obj.Spkevent(6,:);
    obj.Spkevent(10,:)=0*obj.Spkevent(10,:)+1;
    obj.Spkevent(8,:)=0*obj.Spkevent(8,:);
    spkwave=cell(2,1);
    cleanspktt=cell(2,1);
    matwspkPCA=cell(2,1);
    
    nbEl=size(obj.elconfig,2);
    
    for el=1:nbEl
        nbPCmax = obj.Params.MCluster.nbPCmax*(obj.elconfig(2,el)-obj.elconfig(1,el)+1);
        nbPCmin = obj.Params.MCluster.nbPCmin*(obj.elconfig(2,el)-obj.elconfig(1,el)+1);

        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        nbCh=length(chfocus);
        %we first identify the "clean" spikes and perform a spatial 
        %whitening of the spike waveforms
        [spkwave{el},cleanspktt{el}]=obj.ClustPreProcessing(chfocus,timeresol,wnbptsinit,wnbpts);
        fprintf(['el #' num2str(el) ': the first ' num2str(size(cleanspktt{el},2)) ' clean spikes will be used for clustering\n']);
        Ttraining=(obj.Spkevent(2,cleanspktt{el}(end))-obj.Spkevent(2,cleanspktt{el}(1)))*10^-6;
        fprintf(['that corresponds roughly to the first ' num2str(round(Ttraining/60)) ' minutes of the recording\n']);
        
        if nbCh>=1
            matwspk=cell(1,nbCh);                                            
            for ch=1:nbCh
                matwspk{ch}=squeeze(spkwave{el}(ch,:,:));
            end
            %the spike waveforms are projected onto the PCs that show a
            %non-gaussian distribution of their coefficient of projection
            matwspkPCA{el}=obj.globalPCA(matwspk,cleanspktt{el},nbPCmin,nbPCmax,chfocus,totalPC);
            fprintf(['el #' num2str(el) ': we''ll use ' num2str(size(matwspkPCA{el},2)) ' PCs as a feature space\n']);
            totalPC=totalPC+size(matwspkPCA{el},2);
        end
    end
    
    if ~isfield(obj.Params.MCluster,'Seed')
        obj.Params.MCluster.Seed
    end
    
    rng(obj.Params.MCluster.Seed);
    maxIte=obj.Params.MCluster.NbIteration;
    IDXglobalbisAll=zeros(size(obj.Spkwave,2),maxIte,'int16');
    ErrorglobalbisAll=zeros(size(obj.Spkwave,2),maxIte,'single');
    
    chfocusEl=cell(1,nbEl);
    elspkttEl=cell(1,nbEl);
    for el=1:nbEl
        chfocusEl{el}=obj.elconfig(1,el):obj.elconfig(2,el);
        elspkttEl{el}=find(ismember(obj.Spkevent(3,:),chfocusEl{el}));
    end
    Ustart=zeros(1,maxIte);
    
    tic
    for el=1:nbEl                    
        chfocus=chfocusEl{el};
        nbCh=length(chfocus);
        elspktt=elspkttEl{el};
        for ite=1:maxIte
            Ustart(ite)=sum(unique(IDXglobalbisAll(:,ite))>0);
        end

        cleanspktt{el}=cleanspktt{el}(1:min(end,100000));
        matwspkPCAElclean=matwspkPCA{el}(cleanspktt{el},:);
        matwspkPCAEl=matwspkPCA{el};
        spkwaveEl=spkwave{el};

        nbmatwpts=size(matwspkPCAElclean,1);
        nbmaxclust = obj.Params.MCluster.nbmaxunits;
        projrange = obj.Params.MCluster.ProjRange;
        nbreplicates = obj.Params.MCluster.KmeansReplicates;%100;
        Fonline = obj.Params.MCluster.FKmeansOnline;
        
        seed0 = obj.Params.MCluster.Seed;
        if ~FKitespace
            nbunitsglobal0 = obj.Params.MCluster.NbKcluster;
            parfor ite2=1:maxIte
                tic
                disp(['iteration #' num2str(ite2)]);
                cleanspkttIte=cleanspktt{el};
                if nbCh>=1
                    nbunitsglobal = nbunitsglobal0;
                    rng(seed0+ite2);
                    [IDXglobal,~]=kmeans(matwspkPCAElclean,nbunitsglobal,'Start','plus','distance','sqEuclidean','replicates',nbreplicates,'emptyaction','singleton','onlinephase',Fonline);
                    
                    nbspk=zeros(1,nbunitsglobal);
                    for u1=1:nbunitsglobal
                        spkevtu1=find(IDXglobal==u1);
                        nbspk(u1)=numel(spkevtu1);
                    end
                    
                    clustID=sort(unique(IDXglobal),'ascend');
                    nbunitsglobal=numel(clustID);
                    for c1=1:nbunitsglobal
                        IDXglobal(IDXglobal==clustID(c1))=c1;
                    end
                    
                    IDXglobaltemp=IDXglobal;
                    IDXglobal=zeros(size(matwspkPCAEl,1),1);
                    IDXglobal(cleanspkttIte)=IDXglobaltemp;
                    
                    %we measure the chi square value of how well the clusters
                    %explained each spike waveform (in the temporal domain) and
                    %assign all spike to their best fitting clusters
                    Errorglobal=getChi2Parallel(IDXglobal,spkwaveEl,cleanspkttIte,projrange);
                    [chi2,IDXglobal]=min(Errorglobal,[],2);
                    
                    
                    IDXglobalbisAll(elspktt,ite2)=Ustart(ite2)+IDXglobal;
                    ErrorglobalbisAll(elspktt,ite2)=chi2;
                end
                disp(['iteration #' num2str(ite2) ': ' num2str(nbunitsglobal) ' K-means clusters ' num2str(toc) ' sec']);
            end
            obj.IDXglobalAll=round(IDXglobalbisAll);
            obj.ErrorglobalAll=single(ErrorglobalbisAll);
        else
            nbunitsglobal0 = 10:10:200;
            parfor ite2 = 1:numel(nbunitsglobal0)
                tic
                disp(['iteration #' num2str(ite2)]);
                cleanspkttIte=cleanspktt{el};
                if nbCh>=1
                    nbunitsglobal = nbunitsglobal0(ite2);
                    rng(seed0);
                    [IDXglobal,~]=kmeans(matwspkPCAElclean,nbunitsglobal,'Start','plus','distance','sqEuclidean','replicates',nbreplicates,'emptyaction','singleton','onlinephase',Fonline);
                    
                    nbspk=zeros(1,nbunitsglobal);
                    for u1=1:nbunitsglobal
                        spkevtu1=find(IDXglobal==u1);
                        nbspk(u1)=numel(spkevtu1);
                    end
                    
                    clustID=sort(unique(IDXglobal),'ascend');
                    nbunitsglobal=numel(clustID);
                    for c1=1:nbunitsglobal
                        IDXglobal(IDXglobal==clustID(c1))=c1;
                    end
                    
                    IDXglobaltemp=IDXglobal;
                    IDXglobal=zeros(size(matwspkPCAEl,1),1);
                    IDXglobal(cleanspkttIte)=IDXglobaltemp;
                    
                    Errorglobal=getChi2Parallel(IDXglobal,spkwaveEl,cleanspkttIte,projrange);
                    [chi2,IDXglobal]=min(Errorglobal,[],2);
                    
                    
                    IDXglobalbisAll(elspktt,ite2)=Ustart(ite2)+IDXglobal;
                    ErrorglobalbisAll(elspktt,ite2)=chi2;
                    [f,x]=ecdf(mean(ErrorglobalbisAll(elspktt,ite2),2));
                    chi2_95(ite2) = x(find(f>=0.95,1,'first'));
                    chi2_95(ite2) = mean(ErrorglobalbisAll(elspktt,ite2));
                end
                disp(['iteration #' num2str(ite2) ': ' num2str(nbunitsglobal) ' K-means clusters ' num2str(toc) ' sec']);
            end
            save([obj.datafolder filesep obj.filename '_chi2vsKite' num2str(el) '.mat'],'chi2_95','nbunitsglobal0');
            figure;
            plot(nbunitsglobal0,chi2_95);
            hold on;
            plot(nbunitsglobal0(2:end),diff(chi2_95),'r');
        end
    end
    toc
    if FKitespace
        save([obj.datafolder filesep obj.filename '_chi2vsKite.mat'],'chi2_95','nbunitsglobal0');
        figure;
        plot(nbunitsglobal0,chi2_95);
        hold on;
        plot(nbunitsglobal0(2:end),diff(chi2_95),'r');
    end
end