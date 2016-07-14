function dist2clust=getChi2Parallel(IDXglobal,spkwave,cleanspkttIte,ProjRange)
%computes the distance of all spikes to all clusters by fitting the cluster
%templates to the spike waveform, adjusting the template amplitude by a
%scaling factor ranging from 1-ProjRange to 1+ProjRange.
    projlowth=1-ProjRange;
    projhighth=1+ProjRange;
    nbunitsglobal=numel(unique(IDXglobal(IDXglobal>0)));
    nbCh=size(spkwave,1);
    nbSpk=size(spkwave,2);
    nbSpkseg=10000;
    dist2clust=zeros(size(spkwave,2),nbunitsglobal,'single');
    clustcenter=zeros(size(spkwave,1)*size(spkwave,3),nbunitsglobal);
    clustNorm=zeros(1,nbunitsglobal);
    for clust=1:nbunitsglobal
        spkevtidxclean=find(IDXglobal(cleanspkttIte)==clust);
        if numel(spkevtidxclean)>1
            clustcenter(:,clust)=reshape(permute(mean(spkwave(:,cleanspkttIte(spkevtidxclean),:),2),[3 1 2]),[nbCh*size(spkwave,3) 1]);
            clustNorm(clust)=sum((clustcenter(:,clust)).^2,1);
        elseif numel(spkevtidxclean)>0
            clustcenter(:,clust)=reshape(permute(spkwave(:,cleanspkttIte(spkevtidxclean),:),[3 1 2]),[nbCh*size(spkwave,3) 1]);
            clustNorm(clust)=sum((clustcenter(:,clust)).^2,1);
        end
    end
    for seg=1:floor(nbSpk/nbSpkseg)+1
        spkidx=((seg-1)*nbSpkseg+1):min(nbSpk,seg*nbSpkseg);
        spkwf=reshape(permute(spkwave(:,spkidx,:),[3 1 2]),[nbCh*size(spkwave(:,spkidx,:),3) size(spkwave(:,spkidx,:),2)])';
        spkwf = single(spkwf);
        nS=size(spkwf,1);
        ndim=size(spkwf,2);
%         if ProjRange ~= 0
            projcoeff=(spkwf*clustcenter)./repmat(clustNorm,[nS 1]);
            projcoeff=max(projcoeff,projlowth);
            projcoeff=min(projcoeff,projhighth);
            for clust=1:nbunitsglobal
    %             spkteml=projcoeff(:,clust)*clustcenter(:,clust)';
                dist2clust(spkidx,clust)=sum((spkwf-projcoeff(:,clust)*clustcenter(:,clust)').^2,2)/ndim;%/clustNorm(clust);%
            end
%         else
%             for clust=1:nbunitsglobal
%                 dist2clust(spkidx,clust)=sum((bsxfun(@minus,spkwf,clustcenter(:,clust)')).^2,2)/ndim;
%             end
%         end
    end    
end