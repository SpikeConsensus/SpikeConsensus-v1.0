function obj=ShowUnitquality(obj,clustSel)
    unitSel=unique(obj.getUnitIDMultiCh(clustSel))';
    delete(get(obj.dialogUnitquality,'Children'));
    nbclust=size(unitSel,2);
    textaxs=axes('Parent',obj.dialogUnitquality,'Visible','off');
    if isfield(obj.Quality,'FalsePosU') 
        if max(unitSel)<=size(obj.Quality.FalsePosU,2)
            FalsePosUAll2=obj.Quality.FalsePosU;%evalin('base','FalsePosUAll2');
            FalseNegUAll2=obj.Quality.FalseNegU;%evalin('base','FalseNegUAll2');
            FalseNegUMatAll2=obj.Quality.FalseNegUMat;%evalin('base','FalseNegUMatAll2');
            FalsePosUMatAll2=obj.Quality.FalsePosUMat;%evalin('base','FalsePosUMatAll2');
            ConfusionUMatAll2=obj.Quality.ConfusionUMat;%evalin('base','ConfusionUMatAll2');
            CCUMatAll2=obj.Quality.CCUMat;
            SNR2=obj.Quality.SNRindexU;
            CCthresh=0.7;
            if nbclust>0
                nbprop=10;
                unitinfo=zeros(nbclust*nbprop,200);
                l=0;
                for u=1:nbclust
                    try
                    l=l+1;
                    str=['cluster #' num2str(obj.getClusterIDMultiCh(unitSel(u))) '/SNR ' num2str(SNR2(unitSel(u)),3)];
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str=['FP+FN: ' num2str((FalsePosUAll2(unitSel(u))+FalseNegUAll2(unitSel(u)))*100,2)];
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str='FP+FN clust# ';
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    FalseSpk=FalseNegUMatAll2+FalsePosUMatAll2;
                    FalseSpk=max(FalseSpk(unitSel(u),:),[],1);
                    maxFalse=sort(unique(FalseSpk(FalseSpk>0 & FalseSpk<1)),'descend');
                    NBclust=[];
                    NBclustfalse=[];
                    for val=1:min(numel(maxFalse),5)
                        [~,clustNB]=find(ismember(FalseSpk(1,:),maxFalse(val)));
                        for cc=1:numel(clustNB)
                            if CCUMatAll2(unitSel(u),clustNB(cc))>CCthresh
                                NBclust=[NBclust clustNB(cc)];
                                NBclustfalse=[NBclustfalse maxFalse(val)];   
                            end
                        end
                    end
                    %add condition on the CC
                    str=[];
                    for cc=1:numel(NBclust)
                        uclust=obj.getClusterIDMultiCh(NBclust(cc));
                        for uc=1:numel(uclust)
                            str=[str num2str(uclust(uc),3) '/'];
                        end
                        str=[str '   '];
                    end
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str=[];
                    for cc=1:numel(NBclustfalse)
                        str=[str num2str(NBclustfalse(cc)*100,2)];
                        for sp=1:(5-size(num2str(NBclustfalse(cc)*100,2),2))
                            str=[str ' '];
                        end
                    end
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str='Pmis clust#: ';
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    FalseSpk=ConfusionUMatAll2;
                    FalseSpk=max(FalseSpk(unitSel(u),:),[],1);
                    maxFalse=sort(unique(FalseSpk(FalseSpk>0 & FalseSpk<1)),'descend');
                    NBclust=[];
                    NBclustfalse=[];
                    for val=1:min(numel(maxFalse),5)
                        [~,clustNB]=find(ismember(FalseSpk(1,:),maxFalse(val)));
                        for cc=1:numel(clustNB)
                            if CCUMatAll2(unitSel(u),clustNB(cc))>CCthresh
                                NBclust=[NBclust clustNB(cc)];
                                NBclustfalse=[NBclustfalse maxFalse(val)]; 
                            end
                        end
                    end
                    str=[];
                    for cc=1:numel(NBclust)
                        uclust=obj.getClusterIDMultiCh(NBclust(cc));
                        for uc=1:numel(uclust)
                            str=[str num2str(uclust(uc),3) '/'];
                        end
                        str=[str '   '];
                    end
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str=[];
                    for cc=1:numel(NBclustfalse)
                        str=[str num2str(NBclustfalse(cc)*100,2)];
                        for sp=1:(5-size(num2str(NBclustfalse(cc)*100,2),2))
                            str=[str ' '];
                        end
                    end
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str=['Pref: ' num2str(obj.Quality.PrefractoryU(unitSel(u))*100,2)];
                    unitinfo(l,1:length(str))=str;
                    l=l+1;
                    str='';
                    unitinfo(l,1:length(str))=str;                                                        
                    catch
                    end
                end
                htext=text('Parent',textaxs,'string',char(unitinfo),'Units','normalized','FontSize',7,'Position',[0.05 0.5],'HorizontalAlignment','left','VerticalAlignment','middle');
            end
        end
    end
end