function [falsespk goodclust]= GroundTruthTest(obj)
    nbmaxunitperch=obj.nbUnitMultiCh;%numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0)));%numel(obj.SpkwaveclustMulti);%size(obj.SpkclustIDMulti,2);
    
    simulunit=1:numel(unique(obj.Spkevent(11,obj.Spkevent(11,:)>0)));
    nbsimulunit=numel(simulunit);
    
    falsespk=zeros(6,nbsimulunit);
    goodclust=zeros(6,nbsimulunit);
    
    idxgroundtruth=find(obj.Spkevent(11,:)>0);% 
    idxoverlap=find(diff(obj.Spkevent(2,idxgroundtruth))<600);
    idxoverlapgroundtruth=idxgroundtruth(unique([idxoverlap idxoverlap+1 find(obj.Spkevent(12,idxgroundtruth)>0)]));
    
    unitnum=0;
    Nclust = 0;
    Uall=[];
    %for ch=1:nbCh
        for unitID = 1:nbmaxunitperch %Usel=1:nbmaxunitperch
            Usel = obj.getClusterIDMultiCh(unitID)';
            if sum(ismember(obj.SpkeventUM(6,:),Usel))>=1 %special condition for Mcdfiles with groundtruth
            %if obj.SpkclustIDMulti(ch,unitID,1)>0
                unitnum=unitnum+1;
                for Usimul=1:nbsimulunit                    
                    %Usel=squeeze(obj.SpkclustIDMulti(ch,unitID,obj.SpkclustIDMulti(ch,unitID,:)>0));
%                     Uselsimul=obj.SimulUclustID{Usimul};
                    Uselsimul=simulunit(Usimul);
                    Uall=unique([Uall Usel']);
%                     spkSimul=obj.SpkeventSimul(2,ismember(obj.SpkeventSimul(6,:),Uselsimul) & obj.Spkevent(10,:)<2 & obj.SpkeventSimul(2,:)<=max(obj.Spkevent(2,:)));                    
%                     spkEstim=obj.Spkevent(2,ismember(obj.Spkevent(6,:),Usel) & obj.Spkevent(10,:)<2);
%                     
                    %spkSimul=obj.SpkeventSimul(2,ismember(obj.SpkeventSimul(6,:),Uselsimul) & obj.SpkeventSimul(2,:)<=max(obj.Spkevent(2,:)));
                    spkSimul=obj.Spkevent(2,obj.Spkevent(11,:)==Uselsimul | obj.Spkevent(min(12,size(obj.Spkevent,1)),:)==Uselsimul);
                    if numel(unique(obj.SpkeventUM(4,:))) == 1
                        spkEstim=obj.SpkeventUM(2,ismember(obj.SpkeventUM(6,:),Usel));%  & obj.SpkeventUM(4,:)==2);
                    else
                        spkEstim=obj.SpkeventUM(2,ismember(obj.SpkeventUM(6,:),Usel) & obj.SpkeventUM(4,:)==2);
                    end
                    %precision=0.01*obj.Params.detect.Tcensoredfactor*obj.dwsamplefactor*size(obj.SpkwaveclustMulti{1},3)/obj.samplingrate*10^6;
                    precision=2*obj.Params.detect.Tcensoredfactor*obj.dwsamplefactor*size(obj.Spkwave,3)/obj.Params.fileinfo.samplingrate*10^6;
                    NspkSimul=size(spkSimul,2);
                    NspkEstim=size(spkEstim,2);
                    for tspk=1:size(spkEstim,2)
                        if ~isempty(spkSimul(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision))
                            spkSimul(find(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision,1,'first'))=[];
                            spkEstim(tspk)=0;
                            if numel(spkSimul(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision))>1
                                keyboard
                            end
                        end
                    end
                    spkEstim(spkEstim==0)=[];
                    
%                     falsespk(1,Usimul,unitnum)=size(spkEstim,2);%/NspkEstim;%false positive
%                     falsespk(2,Usimul,unitnum)=size(spkSimul,2);%/NspkSimul;%false negative
                    falsespk(3,Usimul)=NspkSimul;
                    if NspkEstim-size(spkEstim,2)>size(spkEstim,2)
                    %if size(spkEstim,2)+size(spkSimul,2)<falsespk(1,Usimul)+falsespk(2,Usimul)
                        falsespk(1,Usimul)=falsespk(1,Usimul)+size(spkEstim,2);%/NspkEstim;%false positive
                        falsespk(2,Usimul)=falsespk(2,Usimul)+NspkEstim-size(spkEstim,2);%/NspkSimul;%true positive
                        falsespk(3,Usimul)=NspkSimul;
                        Nclust = Nclust+1;
                        goodclust(Nclust,Usimul)=unitID;%Usel;
                        
                        spkSimul=obj.Spkevent(2,idxoverlapgroundtruth(obj.Spkevent(11,idxoverlapgroundtruth)==Uselsimul | obj.Spkevent(12,idxoverlapgroundtruth)==Uselsimul));
                        spkEstim=obj.SpkeventUM(2,ismember(obj.SpkeventUM(6,:),Usel));
                        %precision=0.01*obj.Params.detect.Tcensoredfactor*obj.dwsamplefactor*size(obj.SpkwaveclustMulti{1},3)/obj.samplingrate*10^6;
                        precision=2*obj.Params.detect.Tcensoredfactor*obj.dwsamplefactor*size(obj.Spkwave,3)/obj.Params.fileinfo.samplingrate*10^6;
                        NspkSimul=size(spkSimul,2);
                        NspkEstim=size(spkEstim,2);
                        for tspk=1:size(spkEstim,2)
                            if ~isempty(spkSimul(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision))
                                spkSimul(find(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision,1,'first'))=[];
                                spkEstim(tspk)=0;
                                if numel(spkSimul(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision))>1
                                    keyboard
                                end
                            end
                        end
                        spkEstim(spkEstim==0)=[];
                        
                        falsespk(5,Usimul)=falsespk(5,Usimul)+NspkEstim-size(spkEstim,2);%/NspkSimul;%true positive
                        falsespk(6,Usimul)=NspkSimul;
                    end
                end
%                 [failure(1,unitnum),failure(3,unitnum)]=min(squeeze(falsespk(1,:,unitnum)));
%                 [failure(2,unitnum),failure(4,unitnum)]=min(squeeze(falsespk(2,:,unitnum)));
                %failure(3,unitnum)=abs(mean(obj.SpkeventUM(5,ismember(obj.SpkeventUM(6,:),Usel))));
            %end
            end
        end
    %end    
%     unitnum=unitnum+1;
%     spkSimul=obj.SpkeventSimul(2,ismember(obj.SpkeventSimul(6,:),[1:101]) & obj.SpkeventSimul(2,:)<=max(obj.Spkevent(2,:)));
%     spkEstim=obj.Spkevent(2,ismember(obj.Spkevent(6,:),[1:101]));
%     precision=0.01*obj.dwsamplefactor*size(obj.SpkwaveclustMulti{1},3)/obj.samplingrate*10^6;
%     NspkSimul=size(spkSimul,2);
%     NspkEstim=size(spkEstim,2);
%     for tspk=1:size(spkEstim,2)
%         if ~isempty(spkSimul(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision))
%             spkSimul(find(spkSimul>=spkEstim(tspk)-precision & spkSimul<=spkEstim(tspk)+precision,1,'first'))=[];
%             spkEstim(tspk)=0;
%         else
%             %keyboard
%         end
%     end
%     spkEstim(spkEstim==0)=[];
% 
%     failure(1,unitnum)=size(spkEstim,2)/NspkEstim;%/NspkEstim;%false positive
%     failure(2,unitnum)=size(spkSimul,2)/NspkSimul;%/NspkSimul;%false negative
% 
%     if unitnum+1<=size(failure,2)
%         failure(:,unitnum+1:end)=[];
%     end
end