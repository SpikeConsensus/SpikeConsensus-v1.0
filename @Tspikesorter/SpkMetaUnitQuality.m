function obj=SpkMetaUnitQuality(obj,Funit)
% measures the probabilities of misclassification between clusters (or
% single units if Funit is true) from the probability of spikes to be
% co-localized in the same K-means clusters across all N oterations saved 
% in IDXgloalAll.

    ConfusionUMat=cell(1,size(obj.elconfig,2));
    FalsePosU=cell(1,size(obj.elconfig,2));
    FalseNegU=cell(1,size(obj.elconfig,2));
    FalsePosUMat=cell(1,size(obj.elconfig,2));
    FalseNegUMat=cell(1,size(obj.elconfig,2));
    CCUMat=cell(1,size(obj.elconfig,2));
    SNRindex=cell(1,size(obj.elconfig,2));
    waveformChMAX=cell(size(obj.elconfig,2),1);
    nbEl=size(obj.elconfig,2);
    if nargin<2
        Funit=false;
    end
    if Funit
        clust2UnitID=zeros(1,numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0))));
        try
        for clust=1:numel(clust2UnitID)
            clust2UnitID(clust)=obj.getUnitIDMultiCh(clust);
        end
        catch
            keyboard
        end
    else
        clust2UnitID=1:numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0)));
    end
    
    nbunits=numel(unique(clust2UnitID));
    Trefractory = 2;%ms
    if Funit
        obj.Quality.PrefractoryU=zeros(1,nbunits);
        for u=1:nbunits
            ufocus=obj.getClusterIDMultiCh(u);
            obj.Quality.PrefractoryU(u)=obj.RefractViolations(Trefractory,ufocus);
        end
    else
        obj.Quality.PrefractoryC=zeros(1,nbunits);
        for u=1:nbunits
            obj.Quality.PrefractoryC(u)=obj.RefractViolations(Trefractory,u);
        end
    end
    
    spkwave=obj.NoiseWhitening();
    
    maxIte = min(size(obj.IDXglobalAll,2),200);

    for el=1:nbEl
        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        cleanspktt=find(ismember(obj.Spkevent(3,:),chfocus) & ismember(obj.Spkevent(10,:),[0 1 2 3 4 5]) & obj.Spkevent(6,:)>=0);
        cleanspktt2=find(ismember(obj.Spkevent(3,:),chfocus) & ismember(obj.Spkevent(10,:),[0 1 2 3 4 5]) & obj.Spkevent(6,:)>=0);
        IDXglobalbisAll=obj.IDXglobalAll(:,1:maxIte);
        nbIterKmean=size(IDXglobalbisAll,2);
        cleanspktt=cleanspktt(1:min(end,size(IDXglobalbisAll,1)));
        IDXglobalcleanbis=IDXglobalbisAll(cleanspktt,:);
        IDXglobalclean=obj.Spkevent(6,cleanspktt);
        IDXglobalclean2=obj.Spkevent(6,cleanspktt2);
        
        IDXglobalclean(IDXglobalclean<0) = max(IDXglobalclean);

        IDXglobalclean(IDXglobalclean>0)=clust2UnitID(IDXglobalclean(IDXglobalclean>0));
        IDXglobalclean(IDXglobalclean>0)=IDXglobalclean(IDXglobalclean>0)-min(IDXglobalclean(IDXglobalclean>0))+1;
        IDXglobalclean2(IDXglobalclean2>0)=clust2UnitID(IDXglobalclean2(IDXglobalclean2>0));
        IDXglobalclean2(IDXglobalclean2>0)=IDXglobalclean2(IDXglobalclean2>0)-min(IDXglobalclean2(IDXglobalclean2>0))+1;

        ufocus=sort(unique(IDXglobalclean2),'ascend');
        nbunitsglobal=numel(ufocus);
        
        if max(ufocus)~=numel(ufocus)
            error(['wrong number of units detected on electrode #' num2str(el)]);
        end

        nbspk=zeros(1,nbunitsglobal);
        for u1=1:nbunitsglobal
            spkevtu1=find(IDXglobalclean==ufocus(u1));
            nbspk(u1)=numel(spkevtu1);
        end
        tic
        
        wnbptsinit=size(spkwave,3);
        wnbpts=floor(wnbptsinit*obj.Params.detect.Tcensoredfactor);
        if wnbpts>floor(wnbptsinit/2)+1
            wnbpts=floor(wnbpts/2);
        end                
        tkidxfit=(floor(wnbptsinit/2)+1-wnbpts):(floor(wnbptsinit/2)+1+wnbpts);
        clust1=cell(1,nbunitsglobal);

        SNRindex{el}=zeros(1,nbunitsglobal);
        for u=1:nbunitsglobal
            idx=find(IDXglobalclean2==ufocus(u));
            if numel(idx)>1
                clust1{u}=mean(squeeze(spkwave(chfocus,cleanspktt2(idx),tkidxfit)),2);
                stdclust1=std(squeeze(spkwave(chfocus,cleanspktt2(idx),tkidxfit)),1,2);
            elseif ~isempty(idx)
                clust1{u}=spkwave(chfocus,cleanspktt2(idx),tkidxfit);
                stdclust1=clust1{u};
            end
            if ~isempty(idx)
                SNRindex{el}(u)=(max(max(((clust1{u}.^2)./stdclust1.^2))))^0.5;%(max(max((clust1{u}.^2)./(stdclust1.^2))))^0.5;
            end
        end
        
        waveformChMAX{el}=zeros(nbunitsglobal,size(clust1{1},3));
        for u=1:nbunitsglobal
            spkshape=squeeze(clust1{u});
            [chmax,~]=find(abs(spkshape)==max(max(abs(spkshape))));
            waveformChMAX{el}(u,:)=spkshape(chmax,:);
        end
        
        CCUMat{el}=zeros(nbunitsglobal);
        for u1=1:numel(ufocus)-1
            for u2=u1:numel(ufocus)-1
                CCUMat{el}(u1,u2)=sum(sum(clust1{u1}.*clust1{u2}))/(sum(sum(clust1{u1}.^2))*sum(sum(clust1{u2}.^2)))^0.5;%sum(sum(clust1{u1}.*clust1{u2}))/max(sum(sum(clust1{u1}.^2)),sum(sum(clust1{u2}.^2)));
                CCUMat{el}(u2,u1)=CCUMat{el}(u1,u2);
            end
        end
        
        

        ConfusionUMat{el}=zeros(nbunitsglobal);
        FalsePosUMat{el}=zeros(nbunitsglobal);
        FalseNegUMat{el}=zeros(nbunitsglobal);
        FalsePosU{el}=zeros(1,nbunitsglobal);
        FalseNegU{el}=zeros(1,nbunitsglobal);
        for i=1:nbIterKmean
            nbKclust=max(IDXglobalcleanbis(:,i));
            for ku=1:nbKclust
                spkcoloc=find(IDXglobalcleanbis(:,i)==ku);
                ucoloc=sort(unique(IDXglobalclean(spkcoloc)),'ascend');
                sumUcoloc=zeros(1,numel(ucoloc));
                for u1=1:numel(ucoloc)
                    sumUcoloc(u1)=sum(IDXglobalclean(spkcoloc)==ucoloc(u1));
                end
                nbucolocU=numel(ucoloc);
                if nbucolocU>1
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
                end
                for u1=1:numel(ucoloc)
                    Fpos=sumUcoloc(u1);
                    Fneg=sum(sumUcoloc)-sumUcoloc(u1);
                    try
                    if Fpos<Fneg
                        FalsePosU{el}(ucoloc(u1))=FalsePosU{el}(ucoloc(u1))+Fpos/(nbspk(ucoloc(u1))*nbIterKmean);
                    elseif Fneg<Fpos
                        FalseNegU{el}(ucoloc(u1))=FalseNegU{el}(ucoloc(u1))+Fneg/(nbspk(ucoloc(u1))*nbIterKmean);
                    end
                    catch
                        keyboard
                    end
                end
            end
        end
        for u1=1:nbunitsglobal
            ConfusionUMat{el}(u1,u1)=1;
        end 
        toc
    end

    nbunitstotal=numel(unique(clust2UnitID));
    ConfusionUMatAll=zeros(nbunitstotal);
    FalsePosUMatAll=zeros(nbunitstotal,nbunitstotal);
    FalseNegUMatAll=zeros(nbunitstotal,nbunitstotal);
    CCUMatAll=zeros(nbunitstotal);
    nbunitsprevel=0;
    for el=1:nbEl                
        ConfusionUMatAll(nbunitsprevel+1:nbunitsprevel+size(ConfusionUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(ConfusionUMat{el},1))=ConfusionUMat{el};
        FalsePosUMatAll(nbunitsprevel+1:nbunitsprevel+size(FalsePosUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(FalsePosUMat{el},1))=FalsePosUMat{el};
        FalseNegUMatAll(nbunitsprevel+1:nbunitsprevel+size(FalseNegUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(FalseNegUMat{el},1))=FalseNegUMat{el};
        CCUMatAll(nbunitsprevel+1:nbunitsprevel+size(ConfusionUMat{el},1),nbunitsprevel+1:nbunitsprevel+size(ConfusionUMat{el},1))=CCUMat{el};
        nbunitsprevel=nbunitsprevel + size(ConfusionUMat{el},1);
    end
    if Funit
        FalsePosUAll=cell2mat(FalsePosU);
        FalseNegUAll=cell2mat(FalseNegU);
        obj.Quality.ConfusionUMat=ConfusionUMatAll;
        obj.Quality.FalsePosUMat=FalsePosUMatAll;
        obj.Quality.FalseNegUMat=FalseNegUMatAll;
        obj.Quality.FalsePosU=FalsePosUAll;
        obj.Quality.FalseNegU=FalseNegUAll;
        obj.Quality.CCUMat=CCUMatAll;
        SNRindexUAll=cell2mat(SNRindex);
        obj.Quality.SNRindexU=SNRindexUAll;
        waveformChMAXAll=cell2mat(waveformChMAX);
        obj.Quality.SpkshapeMaxU=waveformChMAXAll;
    else
        FalsePosUAll=cell2mat(FalsePosU);
        FalseNegUAll=cell2mat(FalseNegU);
        obj.Quality.ConfusionCMat=ConfusionUMatAll;
        obj.Quality.FalsePosCMat=FalsePosUMatAll;
        obj.Quality.FalseNegCMat=FalseNegUMatAll;
        obj.Quality.FalsePosC=FalsePosUAll;
        obj.Quality.FalseNegC=FalseNegUAll;
        obj.Quality.CCCMat=CCUMatAll;
        SNRindexUAll=cell2mat(SNRindex);
        obj.Quality.SNRindexC=SNRindexUAll;
        waveformChMAXAll=cell2mat(waveformChMAX);
        obj.Quality.SpkshapeMaxU=waveformChMAXAll;
    end
end