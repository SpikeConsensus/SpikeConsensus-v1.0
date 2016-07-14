function matwspk=globalPCA(obj,matwspkch,cleanspktt,nbPCmin,nbPCmax,chfocus,prevPC)
% projects the concatenated spike waveforms on the principal components 
% that show a non gaussian distribution of their coefficients (lilliefors 
% test). The PC features are stored in obj.SpkFeatures. globalPCA is called
% by SpkClusterMultiCh before performing the K-means clusterings.
    nbCh=size(matwspkch,2);    
    matwspk=single(cell2mat(matwspkch));%concatenate spike waveforms across channels        
    wnbpts=floor(size(matwspk,2)/nbCh);
    
    if size(matwspk,1)>=0.1*size(matwspk,2)
        [pcwspk,projpcwspk,~] = princomp(matwspk(cleanspktt,:));        
        lillie=zeros(1,size(projpcwspk,2));
        for kPC=1:min(size(projpcwspk,2),nbPCmax)
            pck=projpcwspk(:,kPC);
            lillie(kPC)=lillietest(pck,0.01);            
        end
        nPC=min(sum(lillie),nbPCmax);
        
        totalPC=0;
        if nPC>0
            PCsel=find(lillie==1);
            PCsel=PCsel(1:min(nbPCmax,end));
            if nPC<nbPCmin
                PCsel=1:nbPCmin;
                nPC=nbPCmin;
            end
            matwspk=matwspk*pcwspk(:,(PCsel));
            for k=1:nPC
                totalPC=totalPC+1;
                obj.SpkFeatures(prevPC+totalPC,chfocus,:)=reshape(pcwspk(:,PCsel(k)),wnbpts,nbCh)';                
            end
        else
            matwspk=[];
        end
    else
        fprintf('no PCA\n')
    end
end