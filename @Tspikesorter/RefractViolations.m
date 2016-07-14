function Prefractory=RefractViolations(obj,Trefractory,unitsel)
    spktimeU=find(ismember(obj.SpkeventUM(6,:),unitsel));
    if ~isempty(spktimeU)     
        fileID=unique(obj.SpkeventUM(1,:));
        spktime=(obj.SpkeventUM(2,obj.SpkeventUM(1,:)==fileID(1) & ismember(obj.SpkeventUM(6,:),unitsel))-obj.SpkeventUM(2,1))*10^-3; 
        offset=0;
        if size(fileID,2)>1
            for k=2:size(fileID,2)
                if ~isempty(spktime)
                    offset=offset+obj.SpkeventUM(2,find(obj.SpkeventUM(1,:)==fileID(k-1),1,'last'))*10^-3;
                end 
                spktime=[spktime offset+(obj.SpkeventUM(2,obj.SpkeventUM(1,:)==fileID(k) & ismember(obj.SpkeventUM(6,:),unitsel))-obj.SpkeventUM(2,1))*10^-3];                            
            end
        end            
        spkISI=diff(spktime);
        Nrefrac=sum(diff(spktime)<Trefractory);

        Prefractory=Nrefrac/numel(spkISI);
    else
        Prefractory=0;
    end
end