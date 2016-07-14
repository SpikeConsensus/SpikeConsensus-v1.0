function obj=BuildNoiseCovMatrix(obj)
% Computes the spatial covariance matrix for each probe based on the noise
% waveforms stored in obj.Noisewave. The covariance matrices are saved in 
% obj.WhiteningMatAllSpatial and used later to do a spatial whitening of
% the detected spike waveforms before clustering.
    obj.Noisewave(:,max(squeeze(obj.Noisewave(1,:,:))')==0,:)=[];    
    obj.Noisewave(isnan(obj.Noisewave)) = 0;
    for el=1:size(obj.elconfig,2)
        nbevtwave=size(obj.Noisewave,2);
        nbch=obj.elconfig(2,el)-obj.elconfig(1,el)+1;
        nbwavpts=size(obj.Noisewave,3);
        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        NoiseCovMat=zeros(nbch,nbch);
        for ch1=1:nbch
            for ch2=ch1:nbch
                noisewavech1=squeeze(double(obj.Noisewave(chfocus(ch1),:,:)/obj.Params.detect.SpkADgain));
                noisewavech2=squeeze(double(obj.Noisewave(chfocus(ch2),:,:)/obj.Params.detect.SpkADgain));
                NoiseCovMat(ch1,ch2)=sum(sum(noisewavech1.*noisewavech2))/(nbwavpts*nbevtwave);
                NoiseCovMat(ch2,ch1)=NoiseCovMat(ch1,ch2);
            end
        end
        obj.WhiteningMatAllSpatial{el}=(NoiseCovMat\eye(size(NoiseCovMat)))^0.5;
    end
end