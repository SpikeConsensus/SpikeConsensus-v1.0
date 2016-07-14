function obj=LoadSortedSpikes(obj,fnickname,Fask)
% loads the sorted spike information (Spkevent, SpkeventUM),
% the K-means cluster labels at every iteration (IDXglobalAll), the
% corresponding chi square values (ErrorglobalAll), the quality measures 
% (Quality) and the correspondance between cluster labels and single unit
% (SpkclustIDMulti)
    if nargin<3
        Fask=false;
    end
    if nargin<2
        fnickname='';
        Fask=false;
    end
    if Fask
        if ischar(obj.datafolder)
            [filname pathfolder]=uigetfile(obj.datafolder);
        else
            [filname pathfolder]=uigetfile(['D:' filesep]);
        end
        fname=[pathfolder filname];
    else
        if isempty(fnickname)
            fname=[obj.datafolder filesep obj.filename '_Spkevent.mat'];
        else
            fname=[obj.datafolder filesep obj.filename '_Spkevent_' fnickname '.mat'];
        end
    end

    obj.SpkeventUM=[];
    obj.MunitChID=[];
    if exist(fname)>0
        load(fname,'info');
        obj.Params.detect.SpkADgain=info.detect.SpkADgain;
        obj.Params.MCluster=info.MCluster;
        if ~isfield(info.MCluster,'nbPCmax')
            obj.Params.MCluster.nbPCmax = obj.Params.MCluster.nbPC;
            obj.Params.MCluster.nbPCmin = 1;
            obj.Params.MCluster = rmfield(obj.Params.MCluster,'nbPC');
        end
        if ~isfield(info.MCluster,'NbKcluster')
            obj.Params.MCluster.NbKcluster = 80;
        end
        if ~isfield(info.MCluster,'KmeansReplicates')
            obj.Params.MCluster.KmeansReplicates = 1;
            obj.Params.MCluster.FKmeansOnline = 'off';
            obj.Params.MCluster.CoreClustersMeth = 'Consistency';
            obj.Params.MCluster.cdfChi2threshnoise = obj.Params.MCluster.cdfChi2thresh;
        end
        
        load(fname,'spkevent');
        obj.Spkevent=spkevent;
        load(fname,'spkeventUM');
        obj.SpkeventUM=spkeventUM;
        
        load(fname,'idxglobalAll');
        obj.IDXglobalAll=idxglobalAll;
        load(fname,'errorglobalAll');
        obj.ErrorglobalAll=errorglobalAll;
        
        load(fname,'quality');
        obj.Quality=quality;
        
        load(fname,'spkclustIDmulti');
        obj.SpkclustIDMulti=spkclustIDmulti;
        load(fname,'spkclustPmisMulti');
        if exist('spkclustPmisMulti')
            obj.SpkclustPmisMulti = spkclustPmisMulti;
        else
            obj.SpkclustPmisMulti = 0.15*ones(1,size(obj.SpkclustIDMulti,1));
        end
        load(fname,'spkwaveclustave');
        obj.SpkwaveclustAve=spkwaveclustave;
        load(fname,'munitchID');
        obj.MunitChID=munitchID;
        load(fname,'mclusterchID');
        obj.MclusterChID=mclusterchID;
        load(fname,'spkfeature');
        obj.SpkFeatures=spkfeature;                       
    else
        warning('no spike sorted yet');
    end
end