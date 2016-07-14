function obj=SpkMetaClassifyMultiCh(obj)
% classifies the spikes which haven't been affected to a cluster after the 
% consensus clustering phase performed by SpkMetaClusterMultiCh. The
% classfication is done by template matching, assigning spikes to the
% cluster whose template explains best the spike waveform in the least
% squares sense. If the chi square value of the best fitting cluster is 
% higher than the chi2 threshold, we consider the spike as a putative overlap of
% spikes and fit recursively the cluster templates to the residual (up to 3
% times). If after this overlap detection phase the chi square value of the
% best fit is still higher than threshold, the spike is classified as
% unsorted and considered as noise.
% Here again, if you want more info or if you really need to understand
% this code, please contact me.

    updatetimewin=obj.Params.MCluster.UpdateTime;           
    cdfchi2thresh=obj.Params.MCluster.cdfChi2thresh;
    cdfchi2threshnoise = obj.Params.MCluster.cdfChi2threshnoise;
    projhighth=1+obj.Params.MCluster.ProjRange;
    projlowth=1-obj.Params.MCluster.ProjRange;
    minnbspkclust=200;
    wnbpts=size(obj.Spkwave,3);
    timeresol=obj.Params.detect.Tcensoredfactor*(floor(size(obj.Spkwave,3))*obj.dwsamplefactor-2)/obj.Params.fileinfo.samplingrate*10^6;
    nbunitstotal=numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0)));
    zerosclustID=nbunitstotal+1;

    obj.Spkevent(6,obj.Spkevent(10,:)>=1)=0;
    obj.SpkeventUM=obj.Spkevent;

    spkwave=obj.NoiseWhitening();

    for el=1:size(obj.elconfig,2)
        disp(['processing el#' num2str(el)]);
        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        nbCh=length(chfocus);

        ErrorglobalbisAll=double(obj.ErrorglobalAll);
        maxIte=min(200,size(ErrorglobalbisAll,2)); 
        ErrorglobalbisAll=ErrorglobalbisAll(:,1:maxIte);
        %we standardize the chi2 error of every iteration
        ErrorglobalbisAll=ErrorglobalbisAll./repmat(mean(ErrorglobalbisAll,1),[size(ErrorglobalbisAll,1) 1])*(mean(reshape(ErrorglobalbisAll,[1 numel(ErrorglobalbisAll)])));

        [f,x]=ecdf(mean(ErrorglobalbisAll(ismember(obj.Spkevent(3,:),chfocus),:),2));
        chi2threshidx=find(f>cdfchi2thresh,1,'first');
        chi2thresh=(round(x(chi2threshidx)/0.001)+1)*0.001;
        chi2threshidx=find(f>cdfchi2threshnoise,1,'first');
        chi2threshnoise=(round(x(chi2threshidx)/0.001)+1)*0.001;
        
        disp(['chi2 threshold overlap= ' num2str(chi2thresh)]);
        disp(['chi2 threshold noise= ' num2str(chi2threshnoise)]);
        
        ErrorglobalbisAll=[];

        spktoclass=find(obj.Spkevent(10,:)>=1 & ismember(obj.Spkevent(3,:),chfocus));
        if ~isempty(spktoclass)
            spktoclass=sort(spktoclass,'ascend');
            
            ufocus=[sort(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0 & ismember(obj.Spkevent(3,:),chfocus))),'ascend') zerosclustID];

            nbunits=length(ufocus)-1;
            ufocusChID=zeros(1,nbunitstotal+1);
            ufocusChID(1:nbunits)=obj.MclusterChID(ufocus(1:end-1));
            
            spktt=spktoclass;
            spktt=[spktt;-1*ones(size(spktt))];
            disp([num2str(length(spktt)) 'spikes left to sort']);
            spkcnt=size(spktt,2);
            
            derivT = diff(obj.Spkevent(2,:));
            Tlength = sum(derivT(derivT>0));
            derivT = [];
            nbspkrefresh=floor(sum(obj.Spkevent(4,spktt(1,:))==el)/((Tlength)*10^-6/updatetimewin))+1;

            %we tag the spikes which are detected as overlapping temporally
            %ie with a interspike interval between 1 and 3 ms
            difftt=[diff(obj.Spkevent(2,spktt(1,:))) 3*timeresol];

            %in case of Spiracle file: we avoid that the time of the
            %last spike in one trial interfere with the first spike of
            %the next one (it's very unlikely but one never know...)
            diffFilenum=[diff(obj.Spkevent(1,spktt(1,:))) 0];
            difftt(diffFilenum>0)=3*timeresol;

            idxoverlap=find(difftt(1:end-1)<timeresol);
            idxoverlap=[];
            spktt(2,idxoverlap)=max(1,floor(min(abs(difftt(idxoverlap)),abs(difftt(idxoverlap)))/10^6*obj.Params.fileinfo.samplingrate/obj.dwsamplefactor));
            spktt(2,idxoverlap+1)=max(spktt(2,idxoverlap+1),0);

            spkttgood=spktt(1,spktt(2,:)==-1);
            spkcntgood=size(spkttgood,2);
            spkttoverlap=spktt(1,spktt(2,:)>0);
            spkcntoverlap=size(spkttoverlap,2);


            delta=floor(3*wnbpts/2);%3*floor(wnbpts/2);
            deltacorrect=wnbpts;%floor(wnbpts/2);
            Clusters=cell(2*delta+1,nbunits+1);
            VarClusters=cell(2*delta+1,nbunits+1);
            for tk=1:2*delta+1
                for u=1:nbunits
                    Clusters{tk,u}=double(zeros(nbCh,1,2*delta));
                    VarClusters{tk,u}=double(zeros(nbCh,1,2*delta));                                                   
                end
                Clusters{tk,nbunits+1}=double(zeros(nbCh,1,2*delta));  
                VarClusters{tk,nbunits+1}=double(zeros(nbCh,1,2*delta));
            end
            varesidual=0;
            ressqr=0;
            sumsqr=0;
            varesidualU=cell(nbunits+1,1);
            nbresidual=0;
            spkwaveclustmulti=cell(nbunits+1,1);
            tkidx=1:wnbpts;
            NarrowFactor=1;
            wnbptswin=floor(2*wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor);
            if floor(wnbpts/2)+1-wnbptswin<0
                wnbptswin=floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor);
            end
            tkidxfit=(floor(wnbpts/2)+1-wnbptswin):(floor(wnbpts/2)+1+wnbptswin);
            BGnoiseidx=tkidxfit;%1:2;
            for  u=1:nbunits
                uidx=find(obj.Spkevent(6,:)==ufocus(u) & obj.Spkevent(10,:)==0);
                if numel(uidx)<2
                    uidx=[uidx uidx];
                end
                try
                spkwaveclustmulti{u}=spkwave(:,uidx,:);
                catch
                    keyboard
                end
                spkwaveclustmulti{u}=spkwaveclustmulti{u}(:,1:min(minnbspkclust,end),:);
                varesidualU{u}=zeros(size(spkwaveclustmulti{u},1),1,size(spkwaveclustmulti{u},3));
                for ch=1:nbCh  
                    mclus=mean(squeeze(spkwaveclustmulti{u}(chfocus(ch),:,:)));
                    clus=[zeros(1,deltacorrect) mclus  zeros(1,deltacorrect)];  
                    vclus=sum((squeeze(spkwaveclustmulti{u}(chfocus(ch),:,:))-repmat(mclus,size(spkwaveclustmulti{u},2),1)).^2)/(size(spkwaveclustmulti{u},2)-1);
                    Varclus=[ones(1,deltacorrect) vclus ones(1,deltacorrect)];
                    ressqr=ressqr+sum(sum((squeeze(spkwaveclustmulti{u}(chfocus(ch),:,tkidxfit))-repmat(mclus(tkidxfit),size(spkwaveclustmulti{u},2),1)).^2));
                    sumsqr=sumsqr+sum(sum((squeeze(spkwaveclustmulti{u}(chfocus(ch),:,tkidxfit))).^2));
                    varesidual=varesidual+sum(sum((squeeze(spkwaveclustmulti{u}(chfocus(ch),:,BGnoiseidx))-repmat(mclus(BGnoiseidx),size(spkwaveclustmulti{u},2),1)).^2));
                    nbresidual=nbresidual+numel(BGnoiseidx)*size(spkwaveclustmulti{u},2);
                    varesidualU{u}(chfocus(ch),1,:)=sum((squeeze(spkwaveclustmulti{u}(chfocus(ch),:,:))-repmat(mclus(:)',size(spkwaveclustmulti{u},2),1)).^2)/(size(spkwaveclustmulti{u},2)-1);
                    if ~isnan(clus)
                        if length(clus)==1
                            clus=[zeros(1,deltacorrect) squeeze(spkwaveclustmulti{u}(chfocus(ch),1,:)) zeros(1,deltacorrect-1)];
                            Varclus=[ones(1,deltacorrect) squeeze(spkwaveclustmulti{u}(chfocus(ch),1,:)) ones(1,deltacorrect-1)];
                        end
                        for tk=delta:-1:0
                            Clusters{delta-tk+1,u}(ch,:,:)=[clus((1+tk):end) zeros(1,tk)];
                            VarClusters{delta-tk+1,u}(ch,:,:)=[Varclus((1+tk):end) ones(1,tk)];
                        end
                        for tk=1:delta
                            Clusters{tk+delta+1,u}(ch,:,:)=[zeros(1,tk) clus(1:(end-tk))];  
                            VarClusters{tk+delta+1,u}(ch,:,:)=[ones(1,tk) Varclus(1:(end-tk))]; 
                        end                                                                               
                    end
                end                                        
            end
            totalerr=ressqr/sumsqr;
            varesidual=varesidual/nbresidual;
            varesidualU{nbunits+1}=varesidual*ones(size(varesidualU{nbunits}));
            spkwaveclustmulti{nbunits+1}=0*spkwaveclustmulti{1}(:,1,:);
            %we create a null cluster to fit the waveforms containing
            %only noise
            for tk=1:2*delta
                Clusters{tk,nbunits+1}=0*Clusters{1,nbunits};
            end
            nbNewSpkwave=zeros(1,nbunits);
            nbOldSpkwave=zeros(1,nbunits);
            timeOldSpkwave=zeros(1,nbunits);
            for  u=1:nbunits
                nbOldSpkwave(u)=find(max(abs(squeeze(spkwaveclustmulti{u}(round(obj.MclusterChID(ufocus(u))),:,:))'))>0,1,'last');
                firstspkidx=obj.Spkevent(2,find(obj.Spkevent(6,:)==ufocus(u),nbOldSpkwave(u),'first'));
                timeOldSpkwave(u)=firstspkidx(end)*10^-6;
            end
            timeNewSpkwave=obj.Spkevent(2,1)*10^-6*ones(1,nbunits);                

            if nbunits>1
                IDXchi=zeros(1,spkcnt);
                valchi2=zeros(1,spkcnt);
                sqrsum2=ones(1,spkcnt);
                nboverlapSpk=0;
                for tseg=1:floor(spkcnt/nbspkrefresh)+1
                    %we first fit the spikes which have been detected
                    %as nonoverlapping spikes.
                    disp([num2str(tseg/(floor(spkcnt/nbspkrefresh)+1)*100) '% through']);
                    ttgood=find(spktt(2,(tseg-1)*nbspkrefresh+1:min(spkcnt,tseg*nbspkrefresh))==-1);
                    %ttgood=ttgood+(tseg-1)*nbspkrefresh;
                    %ttgood=((tseg-1)*nbspkrefresh+1):(min(spkcnt,tseg*nbspkrefresh));
                    if ~isempty(ttgood)
                        ttgood=ttgood+(tseg-1)*nbspkrefresh;
                        if tseg>1
                            for u=1:nbunits
                                if timeNewSpkwave(u)>(timeOldSpkwave(u)+updatetimewin)
                                    disp(['updating cluster waveform of unit#' num2str(u)]);
                                    Spkbefore=obj.Spkevent(2,obj.Spkevent(6,:)==ufocus(u) & obj.Spkevent(2,:)*10^-6<timeNewSpkwave(u));
                                    if numel(Spkbefore)>minnbspkclust
                                        Spkbefore=Spkbefore((numel(Spkbefore)-minnbspkclust):end);
                                        spkwaveclustmulti{u}=spkwave(:,ismember(obj.Spkevent(2,:),Spkbefore),:);%double(obj.NoiseWhitening(obj.Spkwave(:,ismember(obj.Spkevent(2,:),Spkbefore),:)));
                                        for ch=1:nbCh
                                            clus=[zeros(1,deltacorrect) mean(squeeze(spkwaveclustmulti{u}(chfocus(ch),:,:))) zeros(1,deltacorrect)];
                                            if length(clus)==1
                                                clus=[zeros(1,deltacorrect) squeeze(squeeze(spkwaveclustmulti{u}(chfocus(ch),1,:))) zeros(1,deltacorrect-1)];                                                
                                            end
                                            try                                            
                                            for tk=delta:-1:0
                                                Clusters{delta-tk+1,u}(ch,:,:)=[clus((1+tk):end) zeros(1,tk)];
                                            end
                                            for tk=1:delta
                                                Clusters{tk+delta+1,u}(ch,:,:)=[zeros(1,tk) clus(1:(end-tk))];
                                            end
                                            catch 
                                                keyboard
                                            end                                            
                                        end 
                                        nbNewSpkwave(u)=0;
                                        nbOldSpkwave(u)=size(spkwaveclustmulti{u},2);
                                        timeOldSpkwave(u)=timeNewSpkwave(u);
                                    end
                                end
                            end
                        end
                         
                        NarrowFactor=1;
                        disp('processing non-overlapping spikes')                            
                        chi2=1000*ones(length(ttgood),nbunits+1);
                        realchi2=1000*ones(length(ttgood),nbunits+1);
                        tkbest=zeros(1,length(ttgood));
                        Tjitter=4;
                        wnbptswin=floor(2*wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor);
                        if floor(wnbpts/2)+1-wnbptswin<0
                            wnbptswin=floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor);
                        end
                        tkidx=(floor(wnbpts/2)+1-wnbptswin):(floor(wnbpts/2)+1+wnbptswin);%1:wnbpts; 
                        for u=1:nbunits+1
                            for tk=delta+1-Tjitter:delta+1+Tjitter
                                olchi2=min(chi2');
                                templateU=Clusters{tk,u}(:,1,tkidx+deltacorrect);
                                projcoeff=(sum(sum((((spkwave(chfocus,spktt(1,ttgood),tkidx).*repmat(templateU,[1 length(ttgood) 1])))),1),3))'/sum(sum(templateU.^2,1),3);
                                projcoeff=max(projcoeff,projlowth);
                                projcoeff=min(projcoeff,projhighth);
                                spkdifftemp=((spkwave(chfocus,spktt(1,ttgood),tkidx)-repmat(templateU,[1 length(ttgood) 1]).*repmat(projcoeff',[nbCh 1 numel(tkidx)])));                             
                                
                                chi2(:,u)=min(chi2(:,u),(sum(sum(((spkdifftemp).^2),1),3))'/(nbCh*size(tkidx,2)));%/varesidual);
                                tkbest(min(chi2')~=olchi2)=tk;
                            end
                        end
                        sqrsum2(ttgood)=(sum(sum((((spkwave(chfocus,spktt(1,ttgood),tkidx))).^2),1),3))'/(nbCh*size(tkidx,2));%/varesidual;  

                        [valchi2(ttgood),ubest]=min(chi2(:,1:nbunits)');

                        IDXchi(ttgood)=ufocus(ubest);
                        
                        for u=1:nbunits
                            utt=find((IDXchi(ttgood)==ufocus(u)) & valchi2(ttgood)<=chi2thresh & tkbest==delta+1);
                            nbNewSpkwave(u)=nbNewSpkwave(u)+numel(utt);
                            if ~isempty(utt)
                                timeNewSpkwave(u)=obj.Spkevent(2,spktt(1,ttgood(utt(end))))*10^-6;
                            end
                        end
                        
                        %in case the chi2 is above the threshold, we
                        %try to fit the waveforms by the sum of two
                        %spikes
                        chi2notok=find(valchi2(ttgood)>chi2thresh);
                        eventoverlap=[];
                        if ~isempty(chi2notok)
                            disp(['processing type I overlapping spikes: ' num2str(length(chi2notok)) ' instances'])
                            obj.SpkeventUM(10,spktt(1,ttgood(chi2notok)))=2;
                            obj.Spkevent(10,spktt(1,ttgood(chi2notok)))=2;
                            maxoverlap=3;
                            NarrowFactor=0.8;
                            tkidx=(floor(wnbpts/2)+1-floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)):(floor(wnbpts/2)+1+floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor));
%                             tkidx=(floor(wnbpts/2)+1-floor(2*wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)):(floor(wnbpts/2)+1+floor(2*wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor));
                            spkdiff=spkwave(chfocus,spktt(1,ttgood(chi2notok)),tkidx);
                            spkdifftemp=spkdiff;
                            realchi2=10000*ones(length(chi2notok),nbunits);
                            chi2tkoverlap=zeros(maxoverlap,length(chi2notok));
                            ubestoverlap=zeros(maxoverlap,length(chi2notok));
                            for lap=1:maxoverlap
                                if lap==1
                                    Tshiftoverlap1=3;
                                else
                                    Tshiftoverlap1=floor(wnbpts);
                                end
                                for u=1:nbunits
                                    for tk=delta+1-Tshiftoverlap1:delta+1+Tshiftoverlap1
                                        oldrealchi2=min(realchi2,[],2)';                                                        
                                        templateU=Clusters{tk,u}(:,1,tkidx+deltacorrect);
                                        vartemplateU=VarClusters{tk,u}(:,1,tkidx+deltacorrect);
                                        projcoeff=(sum(sum((((spkdiff.*repmat(templateU,[1 length(chi2notok) 1])))),1),3))'/sum(sum(templateU.^2,1),3);
                                        projcoeff=max(projcoeff,projlowth);
                                        projcoeff=min(projcoeff,projhighth);
                                        spkdiffAll=((spkdiff-repmat(templateU,[1 length(chi2notok) 1]).*repmat(projcoeff',[nbCh 1 numel(tkidx)])));
                                        [realchi2(:,u)]=min(realchi2(:,u),(sum(sum(((spkdiffAll).^2),1),3))'/(nbCh*size(tkidx,2)));%/varesidual);
                                        ubestoverlap(lap,min(realchi2,[],2)'<oldrealchi2)=u;
                                        chi2tkoverlap(lap,min(realchi2,[],2)'<oldrealchi2)=tk-(delta+1);
                                        if sum(min(realchi2,[],2)'<oldrealchi2)>0
                                            spkdifftemp(:,min(realchi2,[],2)'<oldrealchi2,:)=spkdiff(:,min(realchi2,[],2)'<oldrealchi2,:)-repmat(templateU,[1 sum(min(realchi2,[],2)'<oldrealchi2) 1]).*repmat(projcoeff(min(realchi2,[],2)'<oldrealchi2)',[nbCh 1 numel(tkidx)]);
                                        end                                              
                                    end
                                end
                                spkdiff=spkdifftemp;
                                oldrealchi2(oldrealchi2<=chi2thresh)=-(abs(oldrealchi2(oldrealchi2<=chi2thresh)));
                            end
                            realchi2=abs(oldrealchi2);

                            sqrsum2(ttgood(chi2notok))=(sum(sum((((spkwave(chfocus,spktt(1,ttgood(chi2notok)),tkidx))).^2),1),3))'/(nbCh*size(tkidx,2));%/varesidual; 
                            valchi2(ttgood(chi2notok))=realchi2;%chi2';%
                            IDXchi(ttgood(chi2notok))=ufocus(ubestoverlap(1,:));
                            ubestoverlap(ubestoverlap==0)=nbunits+1;

                            eventoverlaptemp=cell(1,maxoverlap);
                            for lap=2:maxoverlap
                                chi2u2=ubestoverlap(lap,:);
                                chi2u1=ubestoverlap(1,:);
                                chi2tk=(chi2tkoverlap(lap,:)-chi2tkoverlap(1,:));
                                chi2ok=chi2notok(valchi2(ttgood(chi2notok))<=chi2thresh & abs((chi2tk))<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1));
                                if ~isempty(chi2ok)
                                    try
                                    eventoverlaptemp{lap}=obj.Spkevent(:,spktt(1,ttgood(chi2ok)));
                                    eventoverlaptemp{lap}(2,:)=obj.Spkevent(2,spktt(1,ttgood(chi2ok)))+chi2tk(realchi2<=chi2thresh & abs(chi2tk)<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1))/(obj.Params.fileinfo.samplingrate/obj.dwsamplefactor)*10^6;
                                    eventoverlaptemp{lap}(6,:)=ufocus(chi2u2(realchi2<=chi2thresh & abs(chi2tk)<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1)));
                                    eventoverlaptemp{lap}(7,:)=ufocusChID(ufocus(chi2u2(realchi2<=chi2thresh  & abs(chi2tk)<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1)))); 
                                    eventoverlaptemp{lap}(8,:)=valchi2(ttgood(chi2ok)); 
                                    eventoverlaptemp{lap}(9,:)=sqrsum2(ttgood(chi2ok));
                                    eventoverlaptemp{lap}(10,:)=4;
                                    obj.SpkeventUM(10,spktt(1,ttgood(chi2ok)))=3;
                                    obj.Spkevent(10,spktt(1,ttgood(chi2ok)))=3;
                                    violationidx=find(chi2tk(realchi2<=chi2thresh & abs(chi2tk)<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1))/(obj.Params.fileinfo.samplingrate/obj.dwsamplefactor)*10^6<2000 & chi2u2(realchi2<=chi2thresh & abs(chi2tk)<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1))==chi2u1(realchi2<=chi2thresh & abs(chi2tk)<(floor(wnbpts*obj.Params.detect.Tcensoredfactor*NarrowFactor)-1)));
                                    obj.SpkeventUM(10,spktt(1,ttgood(chi2ok(violationidx))))=5;
                                    obj.Spkevent(10,spktt(1,ttgood(chi2ok(violationidx))))=5;
                                    eventoverlaptemp{lap}(10,violationidx)=6;
                                    for wf=1:numel(chi2ok)
                                        idxuwf=find(obj.Spkevent(6,:)==eventoverlaptemp{lap}(6,wf));
                                        violidx=find(abs(eventoverlaptemp{lap}(2,wf)-obj.Spkevent(2,idxuwf))<=2000,1,'last');
                                        if ~isempty(violidx)
                                            if obj.Spkevent(10,idxuwf(violidx))<=4
                                                eventoverlaptemp{lap}(6,wf)=nbunits+1;
                                                eventoverlaptemp{lap}(10,wf)=4;
                                            end
                                            %keyboard
                                        end
                                    end
                                    catch
                                        keyboard;
                                    end
                                end
                            end
                            eventoverlap=cell2mat(eventoverlaptemp);
                        end
                        %we finally keep the spike identity of the
                        %spikes for which the chi2 threshold is
                        %satisfied
                        chi2ok=find(valchi2(ttgood)<=chi2threshnoise);
                        obj.SpkeventUM(6,spktt(1,ttgood(chi2ok)))=IDXchi(ttgood(chi2ok));
                        obj.SpkeventUM(7,spktt(1,ttgood(chi2ok)))=ufocusChID(IDXchi(ttgood(chi2ok)));
                        obj.Spkevent(6,spktt(1,ttgood(chi2ok)))=IDXchi(ttgood(chi2ok));
                        obj.Spkevent(7,spktt(1,ttgood(chi2ok)))=ufocusChID(IDXchi(ttgood(chi2ok)));
                        obj.SpkeventUM(8,spktt(1,ttgood(chi2ok)))=valchi2(ttgood(chi2ok));%valchi2SingleFit(ttgood(chi2ok));%
                        obj.SpkeventUM(9,spktt(1,ttgood(chi2ok)))=sqrsum2(ttgood(chi2ok));                            
                        obj.Spkevent(8,spktt(1,ttgood(chi2ok)))=valchi2(ttgood(chi2ok));%valchi2SingleFit(ttgood(chi2ok));
                        obj.Spkevent(9,spktt(1,ttgood(chi2ok)))=sqrsum2(ttgood(chi2ok));                            

                        chi2notok=find(valchi2(ttgood)>chi2threshnoise);
                        obj.SpkeventUM(6,spktt(1,ttgood(chi2notok)))=-1;
                        obj.Spkevent(6,spktt(1,ttgood(chi2notok)))=-1;
                        obj.SpkeventUM(8,spktt(1,ttgood(chi2notok)))=valchi2(ttgood(chi2notok));%valchi2SingleFit(ttgood(chi2notok));%
                        obj.SpkeventUM(9,spktt(1,ttgood(chi2notok)))=sqrsum2(ttgood(chi2notok));
                        obj.Spkevent(8,spktt(1,ttgood(chi2notok)))=valchi2(ttgood(chi2notok));%valchi2SingleFit(ttgood(chi2notok));
                        obj.Spkevent(9,spktt(1,ttgood(chi2notok)))=sqrsum2(ttgood(chi2notok));

                        %and we add the spikes which were hidden 
                        %overlapping spikes 
                        %(= which hadn't been detected as overlapping spikes explicitly)                           
                        if ~isempty(eventoverlap)
                            obj.SpkeventUM=[obj.SpkeventUM eventoverlap];
                            nboverlapSpk=nboverlapSpk+size(eventoverlap,2);
                        end                             
                    end

                    %now we fit the spikes which have been detected in
                    %overlap. This is basically the same code as for
                    %detecting hidden overlap, but we are only
                    %interested here in identifying the spk on which
                    %the window is centered                        
                    ttoverlap=find(spktt(2,(tseg-1)*nbspkrefresh+1:min(spkcnt,tseg*nbspkrefresh))>0);
                    if ~isempty(ttoverlap)
                        warning('type II overlapping spikes detected')                                
                    end                          
                end                
            end
            obj.SpkeventUM(6,obj.SpkeventUM(6,:)==zerosclustID)=-1;
            obj.Spkevent(6,obj.Spkevent(6,:)==zerosclustID)=-1;
        end
    end

    fileID=unique(obj.SpkeventUM(1,:));
    globalidx=[];
    for k=1:size(fileID,2)
        IDdx=find(obj.SpkeventUM(1,:)==fileID(k));
        [~,idxreord]=sort(obj.SpkeventUM(2,IDdx),'ascend');
        globalidx=[globalidx IDdx(idxreord)];                
    end
    obj.SpkeventUM=obj.SpkeventUM(:,globalidx);

    nbtotalunits=numel(unique(obj.SpkeventUM(6,obj.SpkeventUM(6,:)>0)));
    for u=1:nbtotalunits
        idxu=find(obj.SpkeventUM(6,:)==u);
        Uisi=diff(obj.SpkeventUM(2,idxu));
        idxoverlap=find(abs(Uisi)<=(floor(wnbpts*obj.Params.detect.Tcensoredfactor)-1)/(obj.Params.fileinfo.samplingrate/obj.dwsamplefactor)*10^6);
        for tspk=1:numel(idxoverlap)
            if obj.SpkeventUM(10,idxu(idxoverlap(tspk)))<=3 || obj.SpkeventUM(10,idxu(idxoverlap(tspk)+1))<=3
                if obj.SpkeventUM(10,idxu(idxoverlap(tspk)))>=3
                    obj.SpkeventUM(6,idxu(idxoverlap(tspk)))=-1;
                elseif obj.SpkeventUM(10,idxu(idxoverlap(tspk)+1))>=3
                    obj.SpkeventUM(6,idxu(idxoverlap(tspk)+1))=-1;
                else
                    %keyboard;
                end
            end
        end                
    end                   
end                
