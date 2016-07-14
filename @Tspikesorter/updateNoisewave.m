function obj=updateNoisewave(obj,rawdata,winsize,dwfactor)
%update the snippets of noise (ie waveforms that didn't cross a threshold
%defined by 4 times the median of the absolute deviations. When enough
%noise waveforms have been accumulated, the noise covariance matrix is
%computed by obj.BuildNoiseCovMatrix.
    HIGHmean=zeros(1,size(rawdata,1));
    HIGHstd=zeros(1,size(rawdata,1));
    thresh=zeros(1,size(rawdata,1));
    zth_noise = 4;
    for ch=1:size(rawdata,1)
        HIGHmean(ch)=0;
        HIGHstd(ch)=median(abs(obj.RawData(ch,:))/0.6745);
        thresh(ch)=HIGHmean(ch)+zth_noise*HIGHstd(ch);
    end

    nbCh=size(rawdata,1);
    threshold=zeros(size(rawdata));
    for ch=1:nbCh
        threshold(ch,:)=thresh(ch);
    end
    
    %detection of both negative and positive crossings
    cross_vec=diff(sign(abs(rawdata)-threshold),1,2);

    
    [~, spkidx]=find(cross_vec<0);
    if size(cross_vec,1)>1
        spkidx=spkidx';
    end


    spkidxdiff=[diff(spkidx) size(rawdata,2)-spkidx(end)];
    spkidxnoisechunk=find(spkidxdiff(1:end-1)>6*winsize);
    nbwavpts=floor(winsize/dwfactor);
    nbCh=size(rawdata,1);
    nbelectrode=size(obj.Params.fileinfo.electrode,2);
    resolnoisewave=min(floor(winsize/dwfactor));
    Nperdim=100;
    if isempty(obj.Noisewave)
        obj.Noisewave=zeros(nbCh,floor(nbCh/nbelectrode)*resolnoisewave*Nperdim,resolnoisewave,'int16');
        num=0;
    else
        num=find(max(abs(squeeze(obj.Noisewave(1,:,:)))')==0,1,'first')-1;
    end
    maxsize=size(obj.Noisewave,2);

    for k=1:size(spkidxnoisechunk,2)-1
        nbchunk=floor((spkidx(spkidxnoisechunk(k)+1)-spkidx(spkidxnoisechunk(k)))/winsize);
        if nbchunk>5
            for kk=2:nbchunk-3
                idx=spkidx(spkidxnoisechunk(k))+kk*winsize:dwfactor:spkidx(spkidxnoisechunk(k))+(kk+1)*winsize;
                num=num+1;
                if num<=maxsize
                    obj.Noisewave(:,num,:)=int16(obj.RawData(:,idx(1:nbwavpts))*obj.Params.detect.SpkADgain);
                end
            end
        end
    end

    if num>=maxsize
        obj.BuildNoiseCovMatrix;
        obj.Noisewave='done';
    end
end