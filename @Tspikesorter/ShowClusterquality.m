function obj=ShowClusterquality(obj,clustSel)
    delete(get(obj.dialogClustquality,'Children'));
    nbclust=size(clustSel,2);
    textaxs=axes('Parent',obj.dialogClustquality,'Visible','off');
    if isfield(obj.Quality,'FalsePosC')
        FalsePosUAll2=obj.Quality.FalsePosC;%evalin('base','FalsePosCAll2');
        FalseNegUAll2=obj.Quality.FalseNegC;%evalin('base','FalseNegCAll2');
        FalseNegUMatAll2=obj.Quality.FalseNegCMat;%evalin('base','FalseNegCMatAll2');
        FalsePosUMatAll2=obj.Quality.FalsePosCMat;%evalin('base','FalsePosCMatAll2');
        ConfusionUMatAll2=obj.Quality.ConfusionCMat;%evalin('base','ConfusionCMatAll2');
        SNR2=obj.Quality.SNRindexC;
        try
        CCUMatAll2=obj.Quality.CCCMat;
        CCthresh=0.7;        
        if nbclust>0
            nbprop=10;
            unitinfo=zeros(nbclust*nbprop,200);
            l=0;
            for u=1:nbclust
                l=l+1;
                str=['cluster #' num2str(clustSel(u),3) '/SNR ' num2str(SNR2(clustSel(u)),3)];
                unitinfo(l,1:length(str))=str;
                l=l+1;
                str=['FP+FN: ' num2str((FalsePosUAll2(clustSel(u))+FalseNegUAll2(clustSel(u)))*100,2)];
                unitinfo(l,1:length(str))=str;
                l=l+1;
                str='FP+FN clust# ';
                unitinfo(l,1:length(str))=str;
                l=l+1;
                FalseSpk=FalseNegUMatAll2+FalsePosUMatAll2;
                FalseSpk=max(FalseSpk(clustSel(u),:),[],1);
                maxFalse=sort(unique(FalseSpk(FalseSpk>0 & FalseSpk<1)),'descend');
                NBclust=[];
                NBclustfalse=[];
                for val=1:min(numel(maxFalse),5)
                    [~,clustNB]=find(ismember(FalseSpk(1,:),maxFalse(val)));
                    for cc=1:numel(clustNB)
                        if CCUMatAll2(clustSel(u),clustNB(cc))>CCthresh
                            NBclust=[NBclust clustNB(cc)];
                            NBclustfalse=[NBclustfalse maxFalse(val)]; 
                        end
                    end
                end
                %add condition on the CC
                str=[];
                for cc=1:numel(NBclust)
                    str=[str num2str(NBclust(cc),3) '    '];
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
                FalseSpk=max(FalseSpk(clustSel(u),:),[],1);
                maxFalse=sort(unique(FalseSpk(FalseSpk>0 & FalseSpk<1)),'descend');
                NBclust=[];
                NBclustfalse=[];
                for val=1:min(numel(maxFalse),5)
                    [~,clustNB]=find(ismember(FalseSpk(1,:),maxFalse(val)));
                    for cc=1:numel(clustNB)
                        if CCUMatAll2(clustSel(u),clustNB(cc))>CCthresh
                            NBclust=[NBclust clustNB(cc)];
                            NBclustfalse=[NBclustfalse maxFalse(val)]; 
                        end
                    end
                end
                str=[];
                for cc=1:numel(NBclust)
                    str=[str num2str(NBclust(cc),3)  '    '];
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
                str='';
                unitinfo(l,1:length(str))=str; 
            end 
            htext=text('Parent',textaxs,'string',char(unitinfo),'Units','normalized','FontSize',7,'Position',[0.05 0.5],'HorizontalAlignment','left','VerticalAlignment','middle');                               
        end
        catch
            warning('quality measures must be updated');
        end
    end
end