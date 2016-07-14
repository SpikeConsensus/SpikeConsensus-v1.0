function [Params]=SpiracleHeaderInfo(obj)
Params.datafolder=obj.filepath;
Params.ADgain=0.000305;
Params.H0gain=5000;    
if ~exist([obj.filepath filesep obj.filename '.mat'])
else
    load([obj.filepath filesep obj.filename],'Episodstartime');
    if ~exist('Episodstartime','var')
        switch obj.acqsystem
            case 'spiracle16'
            case 'Neuralynx'
                fname=[obj.datafolder filesep 'Events.nev'];
                [Timestamps,TTLs] = Nlx2MatEV(fname, [1 0 1 0 0], 0, 1,[]);
                ttlzero=find(TTLs==0);
                Eptag=Timestamps(ttlzero(1:(end-1))+1);
                Episodstartime=[ones(1,length(Eptag));Eptag];
        end
    end    
end

Params.Epstartime=[unique(Episodstartime(1,:));ones(size(unique(Episodstartime(1,:))))];    
end