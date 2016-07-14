function obj=ReorderClusters(obj)
% reorders the clusters according to their proximity in the dendrogram
% associated to the distance matrix computed from the probability of
% misclassification. Requires the matlab function optimalleaforder to
% define the optimal order of the leaf of the dendrogram.
    
    nbEl=size(obj.elconfig,2);
    outperm=cell(1,2);
    for el=1:nbEl
        chfocus=obj.elconfig(1,el):obj.elconfig(2,el);
        cfocus=sort(find(ismember(obj.MclusterChID,chfocus)),'ascend');
        nbunitsglobal=size(obj.Quality.ConfusionCMat(cfocus,cfocus),1);
        rng('default');
        distUmat=squareform(1-obj.Quality.ConfusionCMat(cfocus,cfocus));
        Z=linkage(distUmat);
        f=figure;
        try
            leafOrder=optimalleaforder(Z,distUmat);
            [H,~,outperm{el}]=dendrogram(Z,nbunitsglobal,'Reorder',leafOrder);                
        catch
            [H,~,outperm{el}]=dendrogram(Z,nbunitsglobal);
        end
        outperm{el}=cfocus(outperm{el});
        delete(f);
    end
    outperm=cell2mat(outperm);
    nbunitsglobal=size(obj.Quality.ConfusionCMat,1);
    newlabel=zeros(1,nbunitsglobal);
    for c=1:nbunitsglobal
        newlabel(c)=find(outperm==c);
    end            

    obj.Quality.ConfusionCMat = obj.Quality.ConfusionCMat(outperm,outperm);
    obj.Quality.FalsePosCMat = obj.Quality.FalsePosCMat(outperm,outperm);
    obj.Quality.FalseNegCMat = obj.Quality.FalseNegCMat(outperm,outperm);
    obj.Quality.FalsePosC = obj.Quality.FalsePosC(outperm);
    obj.Quality.FalseNegC = obj.Quality.FalseNegC(outperm);
    obj.Spkevent(6,obj.Spkevent(6,:)>0) = newlabel(obj.Spkevent(6,obj.Spkevent(6,:)>0));
    obj.SpkeventUM(6,obj.SpkeventUM(6,:)>0) = newlabel(obj.SpkeventUM(6,obj.SpkeventUM(6,:)>0));
    obj.SpkclustIDMulti(obj.SpkclustIDMulti>0) = newlabel(obj.SpkclustIDMulti(obj.SpkclustIDMulti>0));
    obj.SpkclustIDMulti = obj.SpkclustIDMulti(outperm,:);
    obj.MclusterChID = obj.MclusterChID(outperm);
    obj.SpkwaveclustAve = obj.SpkwaveclustAve(outperm);    
    
    % to be sure, we re-compute the quality measures for the single unit
    % instead of reoordering their indexes.
    obj.SpkMetaUnitQuality(true);
end 