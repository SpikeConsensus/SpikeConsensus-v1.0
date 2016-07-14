function obj=SpkMetaClusterMultiCh(obj)
    obj=SpkMetaCluster(obj);
    spkeventtemp = obj.Spkevent(6,:);
    nbclust = numel(unique(spkeventtemp(spkeventtemp>0)));
    obj.Spkevent(6,:) = 0*obj.Spkevent(6,:);
    nbunitstart = 0;
    for clust = 1:nbclust
        spkidx = find(spkeventtemp==clust);
        obj=SpkMetaCluster(obj,spkidx,nbunitstart);
        nbunitstart = numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0)));
    end
end

function obj=SpkMetaCluster(obj,spkidxfocus,nbunitstart) 
% performs the consensus based clustering.
% The basic idea is to perform a hierarchical clustering based on a 
% distance matrix defined as the probability of each pair of spikes to be 
% co-localized in the same K-means clusters across all N iterations of 
% K-means clustering. If you need more explanation or if you really
% want to understand this piece of code, please contact me.
    
    if nargin<2 || isempty(spkidxfocus) || isempty(nbunitstart)
        spkidxfocus = 1:size(obj.Spkevent,2);
        nbunitstart = 0;
    end
    
    Fgreedy=obj.Params.MCluster.Fgreedy;                                 
    nbunitstotal=0;
    goodspkidxAll=[];
    ConfusionUMat=cell(1,size(obj.elconfig,2));
    FalsePosU=cell(1,size(obj.elconfig,2));
    FalseNegU=cell(1,size(obj.elconfig,2));
    FalsePosUMat=cell(1,size(obj.elconfig,2));
    FalseNegUMat=cell(1,size(obj.elconfig,2));
    IDXglobalAll=zeros(size(obj.Spkevent(:,spkidxfocus),2),1);
    ErrorglobalAll=zeros(size(obj.Spkevent(:,spkidxfocus),2),1);
    nbEl=size(obj.elconfig,2);
    nbclust=zeros(1,nbEl);
    spkwaveAll=double(obj.NoiseWhitening());
    for el=1:nbEl
        cdfchi2thresh=obj.Params.MCluster.cdfChi2thresh(min(el,end));
        minUSize=obj.Params.MCluster.minUsize(1);
        Psignith=obj.Params.MCluster.Pmisthreshold(min(el,end));
        disp(['processing el#' num2str(el)]);                
        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        IDXglobalbisAll=obj.IDXglobalAll(spkidxfocus,:);
        ErrorglobalbisAll=double(obj.ErrorglobalAll(spkidxfocus,:));
        maxIte=min(200,size(IDXglobalbisAll,2));
        IDXglobalbisAll=IDXglobalbisAll(:,1:maxIte); 
        ErrorglobalbisAll=ErrorglobalbisAll(:,1:maxIte);
        %we standardize the chi2 error of every iteration
        ErrorglobalbisAll=ErrorglobalbisAll./repmat(mean(ErrorglobalbisAll,1),[size(ErrorglobalbisAll,1) 1])*(mean(reshape(ErrorglobalbisAll,[1 numel(ErrorglobalbisAll)])));
        elspkidx=find(ismember(obj.Spkevent(3,spkidxfocus),chfocus));
        ErrorglobalbisAll=ErrorglobalbisAll(elspkidx,:);
        [f,x]=ecdf(mean(ErrorglobalbisAll,2));
        chi2threshidx=find(f>=cdfchi2thresh,1,'first');

        chi2thresh=x(chi2threshidx);
        disp(['chi2 threshold used for metaclustering = ' num2str(chi2thresh)])        
        IDXglobalbisAll=IDXglobalbisAll(elspkidx,:);
        IDXglobalbisAll(mean(ErrorglobalbisAll,2)>chi2thresh,:)=-1; 
        ErrorglobalbisAll=[];
        
        
        allspkidx=1:size(IDXglobalbisAll,1);
        [notgoodspkidx,~]=find(IDXglobalbisAll<0);
        notgoodspkidx=unique(notgoodspkidx);
        goodspkidx=allspkidx(~ismember(allspkidx,notgoodspkidx));
        
        maxnbgoodspk=100000;%we only use the first 100000. Takes too long otherwise...
        goodspkidx=goodspkidx(1:min(end,maxnbgoodspk));
        nbspktotal=numel(goodspkidx);
        disp([num2str(nbspktotal) ' spikes used for consensus clustering']) 

        IDXglobalbis=double(IDXglobalbisAll(goodspkidx,:));

        obj.SpkclustIDMulti=[];
        obj.MunitChID=[];
        obj.MclusterChID=[];
        obj.SpkclustIDMulti=single(zeros(min(20,obj.Params.MCluster.nbmaxunits),min(20,obj.Params.MCluster.nbmaxunits)));
        obj.MunitChID=single(zeros(1,size(obj.elconfig,2)*obj.Params.MCluster.nbmaxunits)); 
        obj.MclusterChID=single(zeros(1,size(obj.elconfig,2)*obj.Params.MCluster.nbmaxunits)); 

        nbIterKmean=size(IDXglobalbis,2);

        IDXglobal=zeros(size(IDXglobalbis,1),1);        
        normIDXglobalbistemp=(sum(IDXglobalbis.^2,2));
        nbunitsglobal=0;
        remainingspkidx=1:size(IDXglobalbis,1);
        while ~isempty(remainingspkidx)
            possiblematchidx=find(normIDXglobalbistemp==normIDXglobalbistemp(1));
            normcrosscorr=(IDXglobalbis(remainingspkidx(possiblematchidx),:)*IDXglobalbis(remainingspkidx(1),:)')./((normIDXglobalbistemp(possiblematchidx)*normIDXglobalbistemp(1)).^0.5);
            spkevt=possiblematchidx(normcrosscorr==1);
            nbunitsglobal=nbunitsglobal+1;
            IDXglobal(remainingspkidx(spkevt))=nbunitsglobal;
            remainingspkidx(spkevt)=[];
            normIDXglobalbistemp(spkevt)=[];
        end

        nbspk=zeros(1,nbunitsglobal);
        for u1=1:nbunitsglobal
            spkevtu1=find(IDXglobal==u1);
            nbspk(u1)=numel(spkevtu1);
        end                    

        notbigenough=find(nbspk<1);
        cleanspktt=find(~ismember(IDXglobal,notbigenough)); 

        IDXglobalcleanbis=IDXglobalbis(cleanspktt,:);
        IDXglobalclean=IDXglobal(cleanspktt);
        clustID=sort(unique(IDXglobalclean),'ascend');
        nbunitsglobal=numel(clustID);
        for c1=1:nbunitsglobal
            IDXglobalclean(IDXglobalclean==clustID(c1))=c1;
        end

        nbspk=zeros(1,nbunitsglobal);
        for u1=1:nbunitsglobal
            spkevtu1=find(IDXglobalclean==u1);
            nbspk(u1)=numel(spkevtu1);
        end 

        if maxIte<255
            ConfusionUMat{el}=zeros(nbunitsglobal,'uint8');
        else
            ConfusionUMat{el}=zeros(nbunitsglobal,'uint16');
        end

        spkevtu=zeros(1,nbunitsglobal);
        for u1=1:nbunitsglobal
            spkevtu(u1)=find(IDXglobalclean==u1,1,'first');
        end
        IDXglobaltemp=IDXglobalcleanbis(spkevtu,:);
        for kite=1:maxIte
            for ku=1:max(IDXglobaltemp(:,kite))
                spkcoloc=find(IDXglobaltemp(:,kite)==ku);
                ConfusionUMat{el}(spkcoloc,spkcoloc)=ConfusionUMat{el}(spkcoloc,spkcoloc)+1;
            end
        end                
        
        IDXglobalcleanClust=IDXglobalclean;

        nbclustglobal=nbunitsglobal;

        n=size(ConfusionUMat{el},2);
        distUmat = nbIterKmean-ConfusionUMat{el}(tril(true(n),-1));
        distUmat = distUmat(:)';
        Z=linkage(single(distUmat)/nbIterKmean);
        Z=double(Z);
        

        %here we try to estimate the cutoff of the hierarchical
        %clustering that gives the most parsimonious (meta) 
        %clustering of the data
        nbclusttarget=1;%100;%100;
        if nbunitsglobal>nbclusttarget
            clustdepth=2;
            nbclusthierarchy=inf(1,100);
            hcutoffstart=0;
            cutoffstep=0.1;
            hcutoff=0;
            while sum(nbclusthierarchy<=nbclusttarget)==0
                hcutoff=hcutoff+1;
                metaClustID=cluster(Z,'cutoff',hcutoffstart+hcutoff*cutoffstep,'depth',clustdepth);
                nbclusthierarchy(hcutoff)=numel(unique(metaClustID));
            end
            nbclusthierarchy(nbclusthierarchy==inf)=[];
            hcutoffstart=hcutoffstart+(find(nbclusthierarchy>nbclusttarget,1,'last'))*cutoffstep;
            nbclusthierarchy=zeros(1,10);
            cutoffstep=0.01;
            parfor hcutoff=1:10
                metaClustID=cluster(Z,'cutoff',hcutoffstart+(hcutoff-1)*cutoffstep,'depth',clustdepth);
                nbclusthierarchy(hcutoff)=numel(unique(metaClustID));
            end
            hcutoffstart=hcutoffstart+(find(nbclusthierarchy>nbclusttarget,1,'last')-1)*cutoffstep;
            cutoffstep=0.001;
            parfor hcutoff=1:10
                metaClustID=cluster(Z,'cutoff',hcutoffstart+(hcutoff-1)*cutoffstep,'depth',clustdepth);
                nbclusthierarchy(hcutoff)=numel(unique(metaClustID));
            end
            hcutoffstart=hcutoffstart+(find(nbclusthierarchy>nbclusttarget,1,'last')-1)*cutoffstep;

            Hclustcutoff=max(cutoffstep,hcutoffstart);
            metaClustID=cluster(Z,'cutoff',Hclustcutoff,'depth',clustdepth);
            if numel(unique(metaClustID))==1
                metaClustID=[1:numel(metaClustID)]';
            end

            IDXglobalclean=metaClustID(IDXglobalclean);                    
            nbunitsglobal=numel(unique(metaClustID));  

            nbspk=zeros(1,nbunitsglobal);
            for u1=1:nbunitsglobal
                spkevtu1=find(IDXglobalclean==u1);
                nbspk(u1)=numel(spkevtu1);
            end
        else
            metaClustID=1:nbunitsglobal;
        end        


        ConfusionMattemp=ConfusionUMat{el};
        IDXglobalcleantemp=IDXglobalclean;
        metaClustIDtemp=metaClustID;
        nbspktemp=nbspk;
        cleanspktttemp=cleanspktt;                        

        mSmax=numel(obj.Params.MCluster.minUsize);
        disp([num2str(mSmax) ' min size threshold values tested (' num2str(obj.Params.MCluster.minUsize) ')']);
        nbClust=zeros(1,mSmax);
        nbFalse=zeros(1,mSmax);
        Fkeepon=true;

        for mS=1:mSmax+1
            if mS==mSmax+1
                optminSizeidx = find(nbClust == max(nbClust),1,'first');%[~,optminSizeidx]=max(nbClust);
                minUSize=obj.Params.MCluster.minUsize(optminSizeidx);
                obj.Params.MCluster.bestminUsize(el)=minUSize;
                disp(['optimal cluster size threshold = ' num2str(minUSize)]) 
                assignin('base','minCSize',minUSize);
            else
                minUSize=obj.Params.MCluster.minUsize(mS);
            end

            if Fkeepon || mS==mSmax+1
                ConfusionUMat{el}=ConfusionMattemp;
                IDXglobalclean=IDXglobalcleantemp;
                nbspk=nbspktemp;
                cleanspktt=cleanspktttemp;
                metaClustID=metaClustIDtemp;

                notbigenough=find(nbspk<minUSize & nbspk>0);
                bigenough=find(nbspk>=minUSize);
                notbigenoughclustID=find(ismember(metaClustID,notbigenough));
                bigenoughclustID=find(ismember(metaClustID,bigenough));
                if numel(bigenoughclustID)>0
                    [PclosestU,~]=max(ConfusionUMat{el}(notbigenoughclustID,bigenoughclustID),[],2);
                    for c1=1:numel(notbigenoughclustID)
                        AllclosestU=find(ConfusionUMat{el}(notbigenoughclustID(c1),bigenoughclustID)==PclosestU(c1));
                        closestU=AllclosestU(nbspk(metaClustID(bigenoughclustID(AllclosestU)))==min(nbspk(metaClustID(bigenoughclustID(AllclosestU)))));
                        closestU=closestU(1);
                        if double(PclosestU(c1))/nbIterKmean>0.02%Psignith && PclosestU>PclosestUsmall%0.01%
                            IDXglobalclean(IDXglobalcleanClust==notbigenoughclustID(c1))=metaClustID(bigenoughclustID(closestU));
                        else
                            IDXglobalclean(IDXglobalcleanClust==notbigenoughclustID(c1))=0;
                        end
                    end 
                else
                    IDXglobalclean(:) = 1;
                end

                cleanspktt=cleanspktt(IDXglobalclean~=0);
                nbspkremoved=sum(IDXglobalclean==0);                
                IDXglobalcleanbis=IDXglobalbis(cleanspktt,:);
                IDXglobalclean=IDXglobalclean(IDXglobalclean~=0);

                unitID=sort(unique(IDXglobalclean),'ascend');
                nbunitsglobal=numel(unitID);
                for u1=1:nbunitsglobal
                    IDXglobalclean(IDXglobalclean==unitID(u1))=u1;
                end 
                if mS==mSmax+1
                    disp([num2str(nbspkremoved) 'spikes cleaned out'])
                    disp(['metaclustering (step 1): ' num2str(nbunitsglobal) ' units found out of ' num2str(nbclustglobal) ' clusters']); 
                end
            end

            if nbspkremoved>0.2*nbspktotal %nbspkremoved>0.02*nbspktotal
                Fkeepon=false;
            end

            if (nbunitsglobal<=300 && nbspkremoved<=0.2*nbspktotal) || mS == 1 %(nbunitsglobal<=300 && nbspkremoved<=0.02*nbspktotal) || mS == 1
                nbspk=zeros(1,nbunitsglobal);
                for u1=1:nbunitsglobal
                    spkevtu1=find(IDXglobalclean==u1);
                    nbspk(u1)=numel(spkevtu1);
                end


                ConfusionUMat{el}=zeros(nbunitsglobal,'single');
                for i=1:nbIterKmean
                    for ku=1:max(IDXglobalcleanbis(:,i))
                        spkcoloc=find(IDXglobalcleanbis(:,i)==ku);
                        ucoloc=sort(unique(IDXglobalclean(spkcoloc)),'ascend');
                        sumUcoloc=zeros(1,numel(ucoloc));
                        for u1=1:numel(ucoloc)
                            sumUcoloc(u1)=sum(IDXglobalclean(spkcoloc)==ucoloc(u1));
                        end
                        for u1=1:numel(ucoloc)-1
                            for u2=u1+1:numel(ucoloc)
                                ConfusionUMat{el}(ucoloc(u1),ucoloc(u2))=ConfusionUMat{el}(ucoloc(u1),ucoloc(u2))+min(sumUcoloc(u1),sumUcoloc(u2))/((nbspk(ucoloc(u1))+nbspk(ucoloc(u2)))*nbIterKmean);
                                ConfusionUMat{el}(ucoloc(u2),ucoloc(u1))=ConfusionUMat{el}(ucoloc(u1),ucoloc(u2));
                            end
                        end                                
                    end
                end
                for u1=1:nbunitsglobal
                    ConfusionUMat{el}(u1,u1)=1;
                end                                                
                
                IDXglobalCore=0*IDXglobal;
                IDXglobalCore(cleanspktt)=IDXglobalclean;

                %in the next step, either we agglomerate cluster closer
                %than the threshold all at once (if Fgreedy is false)
                %or we proceed step by step, by merging the clusters
                %that are the closest and updating the distance matrix
                %ConfusionUmat at every step. In the first case, we
                %favor continuity of the clusters, in the second, we
                %allow continuity only of the clusters to be merged
                %have a relatively sufficient number of spikes 
                if ~Fgreedy
                    try
                    distUmat=squareform((1-ConfusionUMat{el}));
                    Z=linkage(distUmat);
                    catch
                        keyboard
                    end
                    Z=double(Z);
                    f=figure;
%                     leafOrder=optimalleaforder(Z,distUmat);
                    if nbunitsglobal>1
                        [H,T,outperm]=dendrogram(Z,nbunitsglobal);%,'Reorder',leafOrder);
                    else
                        outperm=1;
                    end
                    delete(f);
                    metaClustID=cluster(Z,'cutoff',1-Psignith,'criterion','distance');
                    outpermmeta=metaClustID(outperm);
                    nbunitglobal=numel(unique(metaClustID));
                    newlabel=zeros(1,nbunitglobal);
                    for mclust=1:nbunitglobal
                        oldlabel=outpermmeta(1);
                        newlabel(oldlabel)=mclust;
                        outpermmeta(outpermmeta==oldlabel)=[];
                    end
                    metaClustID=newlabel(metaClustID);
                    IDXglobalclean=metaClustID(IDXglobalclean);
                else                            
                    FgoOn=true;
                    while FgoOn
                        Psignithtemp=max(max(ConfusionUMat{el}-eye(size(ConfusionUMat{el}))));
                        if Psignithtemp>Psignith
                            nbclustglobal=nbunitsglobal;
                            clustlist=1:nbclustglobal;
                            [clust1,clust2]=find((ConfusionUMat{el}-eye(size(ConfusionUMat{el})))==Psignithtemp);
                            clust2merge=unique([clust1' clust2']);

                            metaClustID=clustlist;
                            metaClustID(max(clust2merge))=min(clust2merge);
                            if max(clust2merge)<nbclustglobal
                                metaClustID(max(clust2merge)+1:nbclustglobal)=(max(clust2merge)+1:nbclustglobal)-1;
                            end
                            IDXglobalclean=metaClustID(IDXglobalclean);
                            nbunitsglobal=numel(unique(metaClustID));
                            if mS==mSmax+1
                                disp(['metaclustering (step 1): ' num2str(nbunitsglobal) ' units found out of ' num2str(nbclustglobal) ' clusters at Psigni = ' num2str(Psignithtemp)]);
                            end

                            
                            nbspk=zeros(1,nbunitsglobal);
                            for u1=1:nbunitsglobal
                                spkevtu1=find(IDXglobalclean==u1);
                                nbspk(u1)=numel(spkevtu1);
                            end
                            ConfusionUMat{el}(max(clust2merge),:)=[];
                            ConfusionUMat{el}(:,max(clust2merge))=[];
                            ConfusionUMat{el}(min(clust2merge),:)=0;
                            ConfusionUMat{el}(:,min(clust2merge))=0;
                            for u1=min(clust2merge):min(clust2merge)
                                for u2=1:nbunitsglobal
                                    for i=1:nbIterKmean
                                        sumu=0;
                                        uC1=unique(IDXglobalcleanbis(IDXglobalclean==u1,i));
                                        uC2=unique(IDXglobalcleanbis(IDXglobalclean==u2,i));
                                        uCcom=uC1(ismember(uC1,uC2));
                                        for Ccom=1:numel(uCcom)
                                            sumu1=sum(IDXglobalcleanbis(IDXglobalclean==u1,i)==uCcom(Ccom));
                                            sumu2=sum(IDXglobalcleanbis(IDXglobalclean==u2,i)==uCcom(Ccom));
                                            sumu=sumu+min(sumu1,sumu2);
                                        end
                                        if u1~=u2
                                            ConfusionUMat{el}(u1,u2)=ConfusionUMat{el}(u1,u2)+sumu/((nbspk(u1)+nbspk(u2)));                                            
                                        else
                                            ConfusionUMat{el}(u1,u2)=ConfusionUMat{el}(u1,u2)+sumu/((nbspk(u1)));
                                        end
                                    end
                                    ConfusionUMat{el}(u1,u2)=ConfusionUMat{el}(u1,u2)/nbIterKmean;
                                    ConfusionUMat{el}(u2,u1)=ConfusionUMat{el}(u1,u2);
                                end
                            end
                            FgoOn=true;
                        else
                            FgoOn=false;
                        end
                    end  
                end

                nbspk=zeros(1,nbunitsglobal);
                for u1=1:nbunitsglobal
                    spkevtu1=find(IDXglobalclean==u1);
                    nbspk(u1)=numel(spkevtu1);
                end

                unitID=sort(unique(IDXglobalclean),'ascend');
                nbunitsglobal=numel(unitID);
                for u1=1:nbunitsglobal
                    IDXglobalclean(IDXglobalclean==unitID(u1))=u1;
                end
                if mS==mSmax+1
                    disp(['metaclustering (step 2): ' num2str(nbunitsglobal) ' units found out of ' num2str(nbclustglobal) ' clusters']);
                end
                
                if mS<mSmax+1
                    nbClust(mS) = nbunitsglobal;                  
                end
            end
        end
        if mSmax>1
            figure;plot(nbClust);
        end


        %one might consider reincluding in here the spikes
        %which had been removed because they were part of
        %too small core clusters and with too small
        %connections to the big ones.

        nbspk=zeros(1,nbunitsglobal);
        for u1=1:nbunitsglobal
            spkevtu1=find(IDXglobalclean==u1);
            nbspk(u1)=numel(spkevtu1);
        end
        nbclust(el)=sum(nbspk>=minUSize);%nbunitsglobal;
        disp([num2str(nbclust(el)) ' units found bigger than minUsize ']);

        ConfusionUMat{el}=zeros(nbunitsglobal,'single');
        FalsePosUMat{el}=zeros(nbunitsglobal,'single');
        FalseNegUMat{el}=zeros(nbunitsglobal,'single');
        FalsePosU{el}=zeros(1,nbunitsglobal,'single');
        FalseNegU{el}=zeros(1,nbunitsglobal,'single');
        for i=1:nbIterKmean
            for ku=1:max(IDXglobalcleanbis(:,i))
                spkcoloc=find(IDXglobalcleanbis(:,i)==ku);
                ucoloc=sort(unique(IDXglobalclean(spkcoloc)),'ascend');
                sumUcoloc=zeros(1,numel(ucoloc));
                for u1=1:numel(ucoloc)
                    sumUcoloc(u1)=sum(IDXglobalclean(spkcoloc)==ucoloc(u1));
                end
                for u1=1:numel(ucoloc)-1
                    for u2=u1+1:numel(ucoloc)
                        ConfusionUMat{el}(ucoloc(u1),ucoloc(u2))=ConfusionUMat{el}(ucoloc(u1),ucoloc(u2))+min(sumUcoloc(u1),sumUcoloc(u2))/((nbspk(ucoloc(u1))+nbspk(ucoloc(u2)))*nbIterKmean);
                        ConfusionUMat{el}(ucoloc(u2),ucoloc(u1))=ConfusionUMat{el}(ucoloc(u1),ucoloc(u2));
                        Fpos=sumUcoloc(u1);
                        Fneg=sumUcoloc(u2);
                        if Fpos<Fneg
                            FalsePosUMat{el}(ucoloc(u1),ucoloc(u2))=FalsePosUMat{el}(ucoloc(u1),ucoloc(u2))+Fpos/(nbspk(ucoloc(u1))*nbIterKmean);
                            FalseNegUMat{el}(ucoloc(u2),ucoloc(u1))=FalseNegUMat{el}(ucoloc(u2),ucoloc(u1))+Fpos/(nbspk(ucoloc(u2))*nbIterKmean);
                        else
                            FalseNegUMat{el}(ucoloc(u1),ucoloc(u2))=FalseNegUMat{el}(ucoloc(u1),ucoloc(u2))+Fneg/(nbspk(ucoloc(u1))*nbIterKmean);
                            FalsePosUMat{el}(ucoloc(u2),ucoloc(u1))=FalsePosUMat{el}(ucoloc(u2),ucoloc(u1))+Fneg/(nbspk(ucoloc(u2))*nbIterKmean);
                        end
                    end
                end
                for u1=1:numel(ucoloc)
                    Fpos=sumUcoloc(u1);
                    Fneg=sum(sumUcoloc)-sumUcoloc(u1);
                    if Fpos<Fneg
                        FalsePosU{el}(ucoloc(u1))=FalsePosU{el}(ucoloc(u1))+Fpos/(nbspk(ucoloc(u1))*nbIterKmean);
                    else
                        FalseNegU{el}(ucoloc(u1))=FalseNegU{el}(ucoloc(u1))+Fneg/(nbspk(ucoloc(u1))*nbIterKmean);
                    end
                end
            end
        end
        for u1=1:nbunitsglobal
            ConfusionUMat{el}(u1,u1)=1;
        end
        
        IDXglobal=0*IDXglobal;
        IDXglobal(cleanspktt)=IDXglobalclean;

        CoreClustID=zeros(1,numel(unique(IDXglobalCore(IDXglobalCore~=0))));
        nbCoreClust=numel(unique(IDXglobalCore(IDXglobalCore~=0)));
        for u1=1:nbCoreClust
            CoreClustID(u1)=unique(IDXglobal(IDXglobalCore==u1));
        end

        wnbptsinit=size(spkwaveAll,3);
        NarrowFactor=1;
        wnbpts=floor(2*wnbptsinit*obj.Params.detect.Tcensoredfactor*NarrowFactor);
        if wnbpts>floor(wnbptsinit/2)+1
            wnbpts=floor(wnbpts/2);
        end
        tkidxfit=(floor(wnbptsinit/2)+1-wnbpts):(floor(wnbptsinit/2)+1+wnbpts);
        spkwave=spkwaveAll(chfocus,:,tkidxfit);
        IDXglobalCoretemp=IDXglobalCore;
        IDXglobalCore=zeros(1,numel(elspkidx));
        IDXglobalCore(goodspkidx)=IDXglobalCoretemp;

        Errorchi2=getChi2Parallel(IDXglobalCore,spkwave(:,elspkidx,:),goodspkidx,obj.Params.MCluster.ProjRange);        
        [bestchi2,bestu1]=min(Errorchi2,[],2);
        bestu1(bestchi2>chi2thresh)=0;
        
        IDXglobalCore=bestu1;%0*bestu1;%
        
        bestu1(bestu1>0)=CoreClustID(bestu1(bestu1>0));
        IDXglobal=bestu1;%0*bestu1;%
        
        Errorchi2global=getChi2Parallel(IDXglobal,spkwave(:,elspkidx,:),goodspkidx,obj.Params.MCluster.ProjRange);
        chi2All=zeros(size(IDXglobal));
        for wf=1:numel(chi2All)
            if IDXglobal(wf)>0
                chi2All(wf)=Errorchi2global(wf,IDXglobal(wf));
            end
        end
        IDXglobal(chi2All>chi2thresh)=0;
                           
        IDXglobalAll(elspkidx((IDXglobal>0)))=nbunitstotal+IDXglobal(IDXglobal>0);
        ErrorglobalAll(elspkidx((IDXglobal>0)))=chi2All(IDXglobal>0);
        nbunitstotal=nbunitstotal+nbunitsglobal;
        goodspkidxAll=[goodspkidxAll elspkidx];    
    end
    spkwave=obj.Spkwave;
    nbCh=size(spkwave,1);
    chfocus=1:nbCh;
    clust1=cell(nbCh,nbunitstotal);
    clustpop=cell(nbunitstotal,nbCh);            
    
    for u=1:nbunitstotal
        spkevtidx=find(IDXglobalAll==u);
        for ch=1:nbCh
            if numel(spkevtidx)>1
                clust1{ch,u}=mean(squeeze(spkwave(ch,spkevtidx,:)));
                clustpop{u,ch}=clust1{ch,u};
            else
                clust1{ch,u}=(spkwave(ch,spkevtidx,:));
                clustpop{u,ch}=clust1{ch,u};
            end
        end
    end

    maxpos=zeros(1,nbunitstotal);
    valmax=zeros(1,nbunitstotal);
    for u=1:nbunitstotal
        chmax=1;
        for ch=1:nbCh
            if max(abs(clust1{ch,u}))>valmax(u)
                valmax(u)=max(abs(clust1{ch,u}));
                chmax=chfocus(ch);
            end
        end
        maxpos(u)=chmax;
    end
    chposunit=maxpos;
    ordposidx=1:nbunitstotal;
    %[chposunit,ordposidx]=sort(maxpos,'ascend');
    
    obj.SpkwaveclustAve = cell(1,nbunitstotal);
    for u=1:nbunitstotal
        obj.SpkwaveclustAve{ordposidx(u)} = cell2mat(clust1(:,ordposidx(u)));
    end
    
    ConfusionUMatAll=zeros(nbunitstotal,'single');
    FalsePosUMatAll=zeros(nbunitstotal,'single');
    FalseNegUMatAll=zeros(nbunitstotal,'single');
    nbunitsprevel=0;
    for el=1:nbEl                
        ConfusionUMatAll(nbunitsprevel+1:nbunitsprevel+size(ConfusionUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(ConfusionUMat{el},1))=ConfusionUMat{el};
        FalsePosUMatAll(nbunitsprevel+1:nbunitsprevel+size(FalsePosUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(FalsePosUMat{el},1))=FalsePosUMat{el};
        FalseNegUMatAll(nbunitsprevel+1:nbunitsprevel+size(FalseNegUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(FalseNegUMat{el},1))=FalseNegUMat{el};
        nbunitsprevel=size(ConfusionUMat{el});
    end
    obj.Quality.ConfusionCMat=ConfusionUMatAll(ordposidx,ordposidx);
    obj.Quality.FalsePosCMat=FalsePosUMatAll(ordposidx,ordposidx);
    obj.Quality.FalseNegCMat=FalseNegUMatAll(ordposidx,ordposidx);
    FalsePosUAll=cell2mat(FalsePosU);
    FalseNegUAll=cell2mat(FalseNegU);
    obj.Quality.FalsePosC=FalsePosUAll(ordposidx);
    obj.Quality.FalseNegC=FalseNegUAll(ordposidx);
    
    obj.Quality.Rating=false(nbunitstotal,4);
    
    IDXglobalAll(IDXglobalAll==0)=-nbunitstotal;
    IDXglobalAll=IDXglobalAll+nbunitstotal;
    for c1=1:nbunitstotal
        IDXglobalAll(IDXglobalAll==c1+nbunitstotal)=find(ordposidx==c1);
    end
    if sum(IDXglobalAll>nbunitstotal)~=0
        error('problem when reordering the units')
    end
    
    obj.SpkclustIDMulti=zeros(nbunitstotal,round(nbunitstotal/2));
    for u=1:nbunitstotal
        obj.SpkclustIDMulti(nbunitstart + u,1)=nbunitstart + u;
        obj.MunitChID(nbunitstart + u)=chposunit(u);
        obj.MclusterChID(nbunitstart + u)=chposunit(u);
    end
    totalnbunits=nbunitstart + nbunitstotal;
    obj.MunitChID(totalnbunits+1:end)=[];
    obj.MclusterChID(totalnbunits+1:end)=[];
    
    obj.Spkevent(6,spkidxfocus)=0*obj.Spkevent(6,spkidxfocus);
    obj.Spkevent(6,spkidxfocus)=nbunitstart + IDXglobalAll;
    obj.Spkevent(10,spkidxfocus)=0*obj.Spkevent(10,spkidxfocus);
    obj.Spkevent(10,spkidxfocus(IDXglobalAll==0))=1;
    obj.Spkevent(8,spkidxfocus)=ErrorglobalAll;
    obj.SpkeventUM(:,spkidxfocus)=obj.Spkevent(:,spkidxfocus);               
end