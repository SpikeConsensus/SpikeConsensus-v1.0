function spkwave=NoiseWhitening(obj,SpkWave)
% whitens the spike waveforms with the sptial covariance matrix of the
% noise computed during the detection procedure. After whitening, the spike
% waveforms are also sptailly masked such as the voltage on channels
% located out of the spatial extent of the spike (saved in Spkmask) is
% replaced by a gaussian noise.
    if nargin<2
        nbwavpts=size(obj.Spkwave,3);
        spkwave=zeros(size(obj.Spkwave),'single');        
        for el=1:size(obj.elconfig,2)
            for tt=1:nbwavpts
                spkwave(obj.elconfig(1,el):obj.elconfig(2,el),:,tt)=(squeeze(single(obj.Spkwave(obj.elconfig(1,el):obj.elconfig(2,el),:,tt))'/obj.Params.detect.SpkADgain)*(obj.WhiteningMatAllSpatial{el}))';                
            end
        end
    else
        nbwavpts=size(SpkWave,3);
        spkwave=zeros(size(SpkWave),'single');
        if size(SpkWave,2)>1
            for el=1:size(obj.elconfig,2)
                for tt=1:nbwavpts
                    spkwave(obj.elconfig(1,el):obj.elconfig(2,el),:,tt)=(squeeze(single(SpkWave(obj.elconfig(1,el):obj.elconfig(2,el),:,tt))'/obj.Params.detect.SpkADgain)*(obj.WhiteningMatAllSpatial{el}))';                        
                end
            end
        else
            for el=1:size(obj.elconfig,2)
                spkwave(obj.elconfig(1,el):obj.elconfig(2,el),1,:)=(squeeze(single(SpkWave(obj.elconfig(1,el):obj.elconfig(2,el),1,:)))'/obj.Params.detect.SpkADgain*(obj.WhiteningMatAllSpatial{el}))';                    
            end
        end
    end
    Xval = unique(obj.Params.fileinfo.ChannelPosX);
    Yval = unique(obj.Params.fileinfo.ChannelPosY);
    nbElperunit=min(min(diff(sort(obj.Params.fileinfo.ChannelPosX(obj.Params.fileinfo.ChannelPosY==Yval(1))))),min(diff(sort(obj.Params.fileinfo.ChannelPosY(obj.Params.fileinfo.ChannelPosX==Xval(1))))));
    if isempty(nbElperunit)
        nbElperunit=1;
    end
    spkwinext=2*nbElperunit;
    Poswin=3*nbElperunit;
    nbCh=size(spkwave,1);
    chwin=cell(1,nbCh);
    for ch=1:nbCh
        chwin{ch}=find(abs(obj.Params.fileinfo.ChannelPosY-obj.Params.fileinfo.ChannelPosY(ch))<=Poswin & abs(obj.Params.fileinfo.ChannelPosX-obj.Params.fileinfo.ChannelPosX(ch))<=Poswin);%chfocus;%
    end
    chwinidx=-3:3;
    wnbptsinit=size(obj.Spkwave,3);
    wnbpts=2*2*floor(wnbptsinit*obj.Params.detect.Tcensoredfactor)+1;
    if wnbpts>wnbptsinit
        wnbpts=2*floor(wnbptsinit*obj.Params.detect.Tcensoredfactor)+1;
    end
    nbwf=size(spkwave,2);
    timewin=(floor(wnbptsinit/2)+1-floor(wnbpts/2)):(floor(wnbptsinit/2)+1+floor(wnbpts/2));
    SNRthresh=4*ones(1,numel(chwinidx));
    SNRthresh(floor(numel(chwinidx)/2)+1)=3;
    spkwavemax=max(squeeze(abs(spkwave(:,:,timewin))),[],3);
    rng(1);   
    for ch=1:nbCh
        spkchidx=[];
        if nargin<2
            spkchidx0=find(obj.Params.fileinfo.ChannelPosX(ch)>=obj.Spkmask(1,:)-spkwinext & obj.Params.fileinfo.ChannelPosX(ch)<=obj.Spkmask(2,:)+spkwinext & ...
                           obj.Params.fileinfo.ChannelPosY(ch)>=obj.Spkmask(3,:)-spkwinext & obj.Params.fileinfo.ChannelPosY(ch)<=obj.Spkmask(4,:)+spkwinext);
        else
            spkchidx0=1:size(spkwave,2);
        end
        for chNB=1:numel(chwin{ch})
            if chwin{ch}(chNB)>=1 && chwin{ch}(chNB)<=nbCh
                if chwin{ch}(chNB)==ch
                    SNRthresh=obj.Params.detect.zthresh-2;
                else
                    SNRthresh=obj.Params.detect.zthresh-1;
                end
                [spkchidxchwin]=find(spkwavemax(chwin{ch}(chNB),spkchidx0)>SNRthresh);               
                spkchidx=[spkchidx spkchidx0(spkchidxchwin)];
            end
        end
        
        spkchidx=spkchidx0;
        
        Sigspkchidx=unique(spkchidx);
        noiseidx=find(~ismember(1:nbwf,Sigspkchidx));
        spkwave(ch,noiseidx,:)=randn([1 numel(noiseidx) wnbptsinit]);
    end    
end