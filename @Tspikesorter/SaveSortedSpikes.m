function obj=SaveSortedSpikes(obj,fnickname,fcomment)
% saves under the same folder as the raw data (filename_Spkevent_fnickname.mat) 
% a mat file containing the sorted spike information (Spkevent, SpkeventUM),
% the K-means cluster labels at every iteration (IDXglobalAll), the
% corresponding chi square values (ErrorglobalAll), the quality measures 
% (Quality) and the correspondance between cluster labels and single unit
% (SpkclustIDMulti)
    if nargin<2
        fnickname='';
        fcomment='';
    end
    if nargin<3
        fcomment='';
    end
    if ~isfield(obj.Params.detect,'filepathgroup')
        obj.Params.detect.filepathgroup={obj.filepath};
    end
    for f = 1:numel(obj.Params.detect.filepathgroup)
        datafolder = getdatafolder(obj,obj.Params.detect.filepathgroup{f});
        if ~exist(datafolder)
            datafolder = getdatafolder(obj,obj.filepath);
            warning('didn''t find the path of the appended files so we''ll save them together in the folder of the current file');
        end
        filename = obj.Params.detect.filepathgroup{f}(find(obj.Params.detect.filepathgroup{f}=='/' | obj.Params.detect.filepathgroup{f}=='\',1,'last')+1:end);
        if isempty(fnickname)
            fname=[datafolder filesep filename '_Spkevent.mat'];
        else
            fname=[datafolder filesep filename '_Spkevent_' fnickname '.mat'];
        end

        info=obj.Params;
        spkevent=obj.Spkevent;
        spkeventUM=obj.SpkeventUM;
        idxglobalAll=obj.IDXglobalAll;
        errorglobalAll=obj.ErrorglobalAll;
        quality=obj.Quality;
        spkclustIDmulti=obj.SpkclustIDMulti;
        spkclustPmisMulti = obj.SpkclustPmisMulti;
        spkfeature=obj.SpkFeatures;

        munitchID=obj.MunitChID;
        mclusterchID=obj.MclusterChID;
        nbunitmultich=obj.nbUnitMultiCh;
        spkwaveclustave=obj.SpkwaveclustAve;

        save(fname,'spkevent','spkeventUM','idxglobalAll','errorglobalAll','quality','spkclustIDmulti','spkclustPmisMulti','spkwaveclustave','munitchID','spkfeature','mclusterchID','nbunitmultich','info','fcomment','-v7.3');
    end
end