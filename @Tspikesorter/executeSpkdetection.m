function obj=executeSpkdetection(obj,num)
% detects negative crossings with a threshold defined by obj.Params.detect.zthresh
% times the median of the absolute deviation. Also updates the noise
% covariance matrix by calling the updateNoisewave method.

    if obj.Params.detect.FGroundtruth
        Vmthresh=-20;
        obj.Vmintra(obj.Vmintra==Vmthresh)=Vmthresh-0.01;
        cross_intra=diff(sign(obj.Vmintra-Vmthresh),1,2);
        
        [cellIntra, crossindIntra]=find(cross_intra>0);
    end
    
    HIGHmean=zeros(1,size(obj.RawData,1));
    HIGHstd=zeros(1,size(obj.RawData,1));
    thresh=zeros(1,size(obj.RawData,1));
    %threshold computed from the median of the absolute deviations 
    for ch=1:size(obj.RawData,1)
        HIGHmean(ch)=0;
        HIGHstd(ch)=median(abs(obj.RawData(ch,:))/0.6745);
        thresh(ch)=HIGHmean(ch)+obj.Params.detect.zthresh*HIGHstd(ch);
    end
    
    %censored period in milliseconds
    timeresol=obj.Params.detect.Tcensoredfactor*(obj.Params.detect.wtime2-obj.Params.detect.wtime1);

    nbdatapts=size(obj.RawData,2);
    nbCh=size(obj.RawData,1);
    threshold=zeros(size(obj.RawData));
    for ch=1:nbCh
        threshold(ch,:)=thresh(ch);
    end
    
    %detection of negative crossings
    cross_vec=diff(sign(-obj.RawData-threshold),1,2);
    [chind, crossind]=find(cross_vec>0);        
    crossind = crossind(:)';
    chind = chind(:)';
    
    %censored period in samples
    resoldiff=floor(timeresol/(1000/obj.Params.fileinfo.samplingrate));
    winmax=resoldiff;

    %spike times are defined as the time of the largest trough when the
    %voltage is below the thresold (or in the resoldiff ms follwoing the
    %threshold crossing)
    if ~isempty(crossind)>0
        crossmax=zeros(size(crossind));
        crossenergy=zeros(size(crossind));
        nbevts=size(crossind,2);
        endtime=nbdatapts;
        for evtime=1:nbevts
            if crossind(evtime)+winmax<=endtime
                winspk=find(diff(sign(obj.RawData(chind(evtime),crossind(evtime)+1:crossind(evtime)+winmax)))~=0,1,'first');
                if isempty(winspk)
                    winspk=winmax;
                end
                th=sign(obj.RawData(chind(evtime),crossind(evtime)+1))*thresh(chind(evtime));
                maxpos=[0 diff(sign(diff(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk)))) 0];
                maxpos2=[0 diff(sign(diff(obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):crossind(evtime)+winspk)))) 0];
                if th>0
                    keyboard
                end
                if max(abs(maxpos))==0
                    maxpos(:)=1;
                end
                if th>0
                    iwin=(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk)>th).*(maxpos~=0);
                    iwin2=(obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):crossind(evtime)+winspk)>th).*(maxpos2~=0);
                else
                    iwin=(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk)<th).*(maxpos~=0);
                    iwin2=(obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):crossind(evtime)+winspk)<th).*(maxpos2~=0);
                end
                [smax,ind]=max(abs(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk).*iwin));
                energy=sum((obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):crossind(evtime)+winspk).*iwin2).^2);
            else
                winspk=find(diff(sign(obj.RawData(chind(evtime),crossind(evtime)+1:endtime)))~=0,1,'first');
                if isempty(winspk)
                    winspk=endtime-crossind(evtime);
                end
                th=sign(obj.RawData(chind(evtime),crossind(evtime)+1))*thresh(chind(evtime));
                maxpos=[0 diff(sign(diff(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk)))) 0];
                maxpos2=[0 diff(sign(diff(obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):crossind(evtime)+winspk)))) 0];
                if th>0
                    keyboard
                end
                if max(abs(maxpos))==0
                    maxpos(:)=1;
                end
                if th>0
                    iwin=(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk)>th).*(maxpos~=0);
                    iwin2=(obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):max(crossind(evtime),1)+winspk)>th).*(maxpos2~=0);
                else
                    iwin=(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk)<th).*(maxpos~=0);
                    iwin2=(obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):max(crossind(evtime),1)+winspk)<th).*(maxpos2~=0);
                end
                [smax,ind]=max(abs(obj.RawData(chind(evtime),crossind(evtime):crossind(evtime)+winspk).*iwin));
                energy=sum((obj.RawData(chind(evtime),max(crossind(evtime)-winspk,1):max(crossind(evtime),1)+winspk).*iwin2).^2);
            end
            crossind(evtime)=crossind(evtime)+ind-1;
            crossmax(evtime)=smax/HIGHstd(chind(evtime));
            crossenergy(evtime)=energy*sign(obj.RawData(chind(evtime),crossind(evtime)));
        end

        [crossind,idxreord]=sort(crossind,'ascend');
        chind=chind(idxreord);
        crossmax=crossmax(idxreord);
        crossenergy=crossenergy(idxreord);

        %we identify the successive crossing detections on each channel which
        %are closer than the censored period (ie resoldiff samples) and       
        %keep the spike times with the largest deflection
        for ch=1:nbCh
            evtime=1;
            schind=find(chind==ch);
            if size(schind,2)>1
                while evtime<size(schind,2)
                    while crossind(schind(evtime+1))-crossind(schind(evtime))<resoldiff && evtime<size(schind,2)...
                            && sign(crossmax(schind(evtime+1)))==sign(crossmax(schind(evtime)))
                        if chind(schind(evtime+1))==ch && chind(schind(evtime))==ch
                            if abs(crossmax(schind(evtime+1)))>=abs(crossmax(schind(evtime)));
                                crossind(schind(evtime))=crossind(schind(evtime+1));
                                chind(schind(evtime))=chind(schind(evtime+1));
                                crossmax(schind(evtime))=crossmax(schind(evtime+1));
                                crossenergy(schind(evtime))=crossenergy(schind(evtime+1));
                            end
                            crossind(schind(evtime+1))=[];
                            chind(schind(evtime+1))=[];
                            crossmax(schind(evtime+1))=[];
                            crossenergy(schind(evtime+1))=[];
                            schind(evtime+1)=[];

                            if evtime+1>size(schind,2)
                                break
                            else
                                schind(evtime+1:end)=schind(evtime+1:end)-1;
                            end
                        else
                            error('unexpected shift in spike selection')
                        end
                    end
                    evtime=evtime+1;
                end
            end
        end

        [crossind,idxreord]=sort(crossind,'ascend');
        chind=chind(idxreord);
        crossmax=crossmax(idxreord);
        crossenergy=crossenergy(idxreord);


        for ch=1:nbCh
            evtime=1;
            schind=find(chind==ch);
            if size(schind,2)>1
                while evtime<size(schind,2)
                    while crossind(schind(evtime+1))-crossind(schind(evtime))<resoldiff && evtime<size(schind,2)...
                            && (sign(crossmax(schind(evtime+1)))~=sign(crossmax(schind(evtime))) || crossmax(schind(evtime+1))==crossmax(schind(evtime)))
                        if chind(schind(evtime+1))==ch && chind(schind(evtime))==ch
                            if (crossmax(schind(evtime+1))<0 && crossmax(schind(evtime+1))+crossmax(schind(evtime))<15) || abs(crossmax(schind(evtime+1)))>=1.5*abs(crossmax(schind(evtime)));
                                crossind(schind(evtime))=crossind(schind(evtime+1));
                                chind(schind(evtime))=chind(schind(evtime+1));
                                crossmax(schind(evtime))=crossmax(schind(evtime+1));
                                crossenergy(schind(evtime))=crossenergy(schind(evtime+1));
                            end
                            crossind(schind(evtime+1))=[];
                            chind(schind(evtime+1))=[];
                            crossmax(schind(evtime+1))=[];
                            crossenergy(schind(evtime+1))=[];
                            schind(evtime+1)=[];

                            if evtime+1>size(schind,2)
                                break
                            else
                                schind(evtime+1:end)=schind(evtime+1:end)-1;
                            end
                        else
                            error('unexpected shift in spike selection')
                        end
                        if crossind(schind(evtime+1))-crossind(schind(evtime))<resoldiff && evtime<size(schind,2)...
                                && (sign(crossmax(schind(evtime+1)))==sign(crossmax(schind(evtime))) || crossmax(schind(evtime+1))==crossmax(schind(evtime)))
                            if chind(schind(evtime+1))==ch && chind(schind(evtime))==ch
                                if abs(crossmax(schind(evtime+1)))>=abs(crossmax(schind(evtime)))
                                    crossind(schind(evtime))=crossind(schind(evtime+1));
                                    chind(schind(evtime))=chind(schind(evtime+1));
                                    crossmax(schind(evtime))=crossmax(schind(evtime+1));
                                    crossenergy(schind(evtime))=crossenergy(schind(evtime+1));
                                end
                                crossind(schind(evtime+1))=[];
                                chind(schind(evtime+1))=[];
                                crossmax(schind(evtime+1))=[];
                                crossenergy(schind(evtime+1))=[];
                                schind(evtime+1)=[];

                                if evtime+1>size(schind,2)
                                    break
                                else
                                    schind(evtime+1:end)=schind(evtime+1:end)-1;
                                end
                            else
                                error('unexpected shift in spike selection')
                            end
                        end
                    end
                    evtime=evtime+1;
                end
            end
        end

        [crossind,idxreord]=sort(crossind,'ascend');
        chind=chind(idxreord);
        crossmax=crossmax(idxreord);
        crossenergy=crossenergy(idxreord);

        for ch=1:nbCh
            evtime=1;
            schind=find(chind==ch);
            if size(schind,2)>1
                while evtime<size(schind,2)
                    while crossind(schind(evtime+1))-crossind(schind(evtime))<resoldiff && evtime<size(schind,2)... %resoldiffoverlap
                            && sign(crossmax(schind(evtime+1)))==sign(crossmax(schind(evtime)))
                        if chind(schind(evtime+1))==ch && chind(schind(evtime))==ch
                            if abs(crossmax(schind(evtime+1)))>=abs(crossmax(schind(evtime)));%abs(crossmax(evtime+1))>=abs(crossmax(evtime))
                                crossind(schind(evtime))=crossind(schind(evtime+1));
                                chind(schind(evtime))=chind(schind(evtime+1));
                                crossmax(schind(evtime))=crossmax(schind(evtime+1));
                                crossenergy(schind(evtime))=crossenergy(schind(evtime+1));
                            end
                            crossind(schind(evtime+1))=[];
                            chind(schind(evtime+1))=[];
                            crossmax(schind(evtime+1))=[];
                            crossenergy(schind(evtime+1))=[];
                            schind(evtime+1)=[];

                            if evtime+1>size(schind,2)
                                break
                            else
                                schind(evtime+1:end)=schind(evtime+1:end)-1;
                            end
                        else
                            error('unexpected shift in spike selection')
                        end
                    end
                    evtime=evtime+1;
                end
            end
        end



        [crossind,idxreord]=sort(crossind,'ascend');
        chind=chind(idxreord);
        crossmax=crossmax(idxreord);
        crossenergy=crossenergy(idxreord);

        filenum=obj.Epstartime(1,num)*ones(size(crossind));
        elind=zeros(size(chind));
        for el=1:size(obj.elconfig,2)
            elind(chind>=obj.elconfig(1,el) & chind<=obj.elconfig(2,el))=el;
        end
                
        %we identify threshold crossing detected on different channels but
        %closer than the censored period (resoldiff). Those that are closer
        %than 3 electrode units are lumped together and we only keep the
        %spike time corresponding to the largest deflection. During this
        %procedure, we keep track of the spatial extent of each detected
        %spike by saving in maskind the X and Y extent of the threshold 
        %crossing which are lumped together
        maskind=zeros(4,numel(crossind));
        maskind(1,:)=obj.Params.fileinfo.ChannelPosX(chind);
        maskind(2,:)=obj.Params.fileinfo.ChannelPosX(chind);
        maskind(3,:)=obj.Params.fileinfo.ChannelPosY(chind);
        maskind(4,:)=obj.Params.fileinfo.ChannelPosY(chind);
        
        fullcleanspktt=[];
        for el=1:size(obj.elconfig,2)
            disp(['Multi Channel Cleaning #' num2str(el)]);
            chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
            nbElperunit=min(min(diff(sort(obj.Params.fileinfo.ChannelPosX(obj.Params.fileinfo.ChannelPosY==obj.Params.fileinfo.ChannelPosY(chfocus(1)))))),min(diff(sort(obj.Params.fileinfo.ChannelPosY(obj.Params.fileinfo.ChannelPosX==obj.Params.fileinfo.ChannelPosX(chfocus(1)))))));
            if isempty(nbElperunit)
                nbElperunit=1;
            end
            nbchel=3*nbElperunit;%length(chfocus);

            cleanspktt=find(ismember(chind,chfocus));
            cleanspktt=sort(cleanspktt,'ascend');
            
            nbspkmerged=1;
            while nbspkmerged>0
                nbspkmerged=0;
                tt=0;
                while tt<size(cleanspktt,2)
                    tt=tt+1;tk=1;
                    while tk<=numel(chfocus) && tt+tk<=size(cleanspktt,2)
                        minchX(1)=maskind(1,cleanspktt(tt));
                        minchX(2)=maskind(1,cleanspktt(tt+tk));
                        maxchX(1)=maskind(2,cleanspktt(tt));
                        maxchX(2)=maskind(2,cleanspktt(tt+tk));
                        minchY(1)=maskind(3,cleanspktt(tt));
                        minchY(2)=maskind(3,cleanspktt(tt+tk));
                        maxchY(1)=maskind(4,cleanspktt(tt));
                        maxchY(2)=maskind(4,cleanspktt(tt+tk));
                        if abs(crossind(cleanspktt(tt+tk))-crossind(cleanspktt(tt)))<resoldiff...%ie <0.25*duration of spk waveforms (1ms for a 4ms time window)
                                && (abs(minchY(1)-maxchY(2))<=nbchel || abs(maxchY(1)-minchY(2))<=nbchel)...
                                && (abs(minchX(1)-maxchX(2))<=nbchel || abs(maxchX(1)-minchX(2))<=nbchel)...
                                && elind(cleanspktt(tt+tk))==elind(cleanspktt(tt))...
                                && filenum(cleanspktt(tt+tk))==filenum(cleanspktt(tt))
                            nbspkmerged=nbspkmerged+1;
                            if abs(crossmax(cleanspktt(tt+tk)))>=abs(crossmax(cleanspktt(tt)))
                                maskind(1,cleanspktt(tt+tk))=min(maskind(1,cleanspktt(tt+tk)),maskind(1,cleanspktt(tt)));
                                maskind(2,cleanspktt(tt+tk))=max(maskind(2,cleanspktt(tt+tk)),maskind(2,cleanspktt(tt)));
                                maskind(3,cleanspktt(tt+tk))=min(maskind(3,cleanspktt(tt+tk)),maskind(3,cleanspktt(tt)));
                                maskind(4,cleanspktt(tt+tk))=max(maskind(4,cleanspktt(tt+tk)),maskind(4,cleanspktt(tt)));
                                cleanspktt(tt)=[];
                                tt=max(1,tt-1);
                                tk=1;
                            else
                                maskind(1,cleanspktt(tt))=min(maskind(1,cleanspktt(tt+tk)),maskind(1,cleanspktt(tt)));
                                maskind(2,cleanspktt(tt))=max(maskind(2,cleanspktt(tt+tk)),maskind(2,cleanspktt(tt)));
                                maskind(3,cleanspktt(tt))=min(maskind(3,cleanspktt(tt+tk)),maskind(3,cleanspktt(tt)));
                                maskind(4,cleanspktt(tt))=max(maskind(4,cleanspktt(tt+tk)),maskind(4,cleanspktt(tt)));
                                cleanspktt(tt+tk)=[];
                                tk=1;
                            end
                        else
                            tk=tk+1;
                        end
                    end
                end
            end
            try
            fullcleanspktt=sort([fullcleanspktt cleanspktt],'ascend');
            catch
                keyboard
            end
        end
        crossind=crossind(fullcleanspktt);
        chind=chind(fullcleanspktt);
        maskind=maskind(:,fullcleanspktt);
        crossmax=crossmax(fullcleanspktt);
        crossenergy=crossenergy(fullcleanspktt);

        [crossind,idxreord]=sort(crossind,'ascend');
        chind=chind(idxreord);
        crossmax=crossmax(idxreord);
        maskind=maskind(:,idxreord);
        crossenergy=crossenergy(idxreord);
        
        nbevts=size(crossind,2);
        spkevtTime=obj.timestamp1+(crossind/obj.Params.fileinfo.samplingrate)*10^6;
        filenum=ones(size(crossind));
        elind=zeros(size(chind));
        unitID=zeros(size(crossind));
        chi2=-1*ones(size(crossind));
        sqrsum=ones(size(crossind));
        
        for el=1:size(obj.elconfig,2)
            elind(chind>=obj.elconfig(1,el) & chind<=obj.elconfig(2,el))=el;
        end
        
        %we fill Spkevent, SpkeventUM and Spkmask with the information we
        %have on each spike
        spkvec=[filenum;spkevtTime;chind;elind;crossmax;unitID;unitID;chi2;sqrsum;zeros(1,numel(spkevtTime));zeros(1,numel(spkevtTime));zeros(1,numel(spkevtTime))];
        obj.Spkevent=[obj.Spkevent spkvec];
        obj.Spkmask=[obj.Spkmask maskind];
        obj.SpkeventUM=obj.Spkevent;
        
        %we save the spike waveforms in Spkwave, aligned on the time of
        %maximal ampliture. Note that after alignement, the waveforms are
        %downsampled to the Nyquist frequency of the high pass filter
        %(obj.dwsamplefactor)
        win=floor(0.5*(obj.Params.detect.wtime2-obj.Params.detect.wtime1)/(1000/obj.Params.fileinfo.samplingrate));
        spkwav=int16(zeros(nbCh,nbevts,floor(2*win/obj.dwsamplefactor)));
        nbwavpts=size(spkwav,3);
        disp(nbevts)
        try
        nbrawpts=size(obj.RawData,2);
        for k=1:nbevts
            idx=crossind(k)-win:obj.dwsamplefactor:crossind(k)+win;
            idx=max(idx,1);idx=min(idx,nbrawpts);
            spkwav(:,k,:)=int16(obj.RawData(:,idx(1:nbwavpts))*obj.Params.detect.SpkADgain);
        end
        catch
            keyboard
        end
        obj.Spkwave=int16(cat(2,obj.Spkwave,spkwav));
        
        %if the noise covariance matrix hasn't been computed yet, we still
        %accumulate snippets of noise waveforms in Noisewave
        if ~ischar(obj.Noisewave)
            obj.updateNoisewave(obj.RawData,2*win,obj.dwsamplefactor);
        end
        
        %if this is a file with ground truth, we label the detected spikes 
        %that are the closest to the intracellular AP as groundtruth spikes
        if obj.Params.detect.FGroundtruth
            APevt=(obj.timestamp1+((crossindIntra)/obj.Params.fileinfo.samplingrate)*10^6);
            if size(APevt,1)>1
                APevt=APevt';
            end
            obj.APevent=[obj.APevent APevt];
            precision=2*obj.Params.detect.Tcensoredfactor*obj.dwsamplefactor*size(obj.Spkwave,3)/obj.Params.fileinfo.samplingrate*10^6;
            for tspk=1:numel(APevt)
                if cellIntra(tspk)>=1
                    idxspk=find(obj.Spkevent(2,:)>=APevt(tspk)-precision & obj.Spkevent(2,:)<=APevt(tspk)+precision);% & obj.Spkevent(4,:)==2);
                    [~,closestspk]=min(abs(obj.Spkevent(2,idxspk)-APevt(tspk)));
                    if obj.Spkevent(11,idxspk(closestspk))==0
                        obj.Spkevent(11,idxspk(closestspk))=cellIntra(tspk);
                    else
                        obj.Spkevent(12,idxspk(closestspk))=cellIntra(tspk);
                    end
                end
            end
            disp([num2str(numel(APevt(cellIntra>=1))-(numel(find(obj.Spkevent(11,:)>=1))+numel(find(obj.Spkevent(12,:)>=1)))) ' intracellular spikes not detected extracellularly']);

        end
    end
end