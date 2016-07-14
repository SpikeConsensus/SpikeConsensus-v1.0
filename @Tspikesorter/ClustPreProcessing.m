function [spkwave,cleanspktt]=ClustPreProcessing(obj,chfocus,timeresol,wnbptsinit,wnbpts)
% selects the "clean" spikes, whitens them with the noise covariance matrix
% and truncate them to a temporal window of wnbpts points centered on the 
% spike time. ClustPreProcessing is called by SpkClusterMultiCh before the
% feature extraction step.

    nbElperunit=min(min(diff(sort(obj.Params.fileinfo.ChannelPosX(obj.Params.fileinfo.ChannelPosY==obj.Params.fileinfo.ChannelPosY(chfocus(1)))))),min(diff(sort(obj.Params.fileinfo.ChannelPosY(obj.Params.fileinfo.ChannelPosX==obj.Params.fileinfo.ChannelPosX(chfocus(1)))))));
    if isempty(nbElperunit)
        nbElperunit=1;
    end
    spkwinext=2*nbElperunit;
    nbCh=numel(chfocus);
    cleanspktt=cell(1,nbCh);
    for ch=1:nbCh
        cleanspktt{ch}=find(obj.Spkevent(3,:)==chfocus(ch));% & abs(obj.Spkevent(5,:))>=6);
    end

    elspktt=find(ismember(obj.Spkevent(3,:),chfocus));
    cleanspktt=1:numel(elspktt);
    disp(length(cleanspktt));
    
    %we won't consider for the initial part of the clustering, spikes that
    %are spatiotemporaly overlapping, ie with masks distant by less than 2 
    %channels and spike times closer than timeresol.
    idxoverlap=[];
    for ch=1:nbCh
        winspktt=find(obj.Params.fileinfo.ChannelPosX(chfocus(ch))>=obj.Spkmask(1,elspktt)-spkwinext & obj.Params.fileinfo.ChannelPosX(chfocus(ch))<=obj.Spkmask(2,elspktt)+spkwinext &...
                      obj.Params.fileinfo.ChannelPosY(chfocus(ch))>=obj.Spkmask(3,elspktt)-spkwinext & obj.Params.fileinfo.ChannelPosY(chfocus(ch))<=obj.Spkmask(4,elspktt)+spkwinext);
        difftt=[diff(obj.Spkevent(2,elspktt(winspktt))) 2*timeresol];
        idxoverlapch=find(difftt(1:end-1)<timeresol);
        idxoverlap=unique([idxoverlap winspktt(idxoverlapch) winspktt(idxoverlapch+1)]);
    end
    cleanspktt(idxoverlap)=[];
    
    %we'll consider only the first obj.Params.MCluster.TrainingMaxnbSpk spikes 
    %at most for the intial Kmeans clustering phase.
    Nspktraining=min(obj.Params.MCluster.TrainingMaxnbSpk,size(obj.Spkevent,2));    
    cleanspktt=cleanspktt(1:min(end,Nspktraining));
    
    %we spatially whiten the spike waveforms and only consider a time
    %window of wnbpts samples centered on the spike time
    spkwave = obj.NoiseWhitening();
    spkwave=spkwave(chfocus,elspktt,(floor(wnbptsinit/2)+1-floor(wnbpts/2)):(floor(wnbptsinit/2)+1+floor(wnbpts/2)));   
end