function obj=LoadSpikeData(obj,fnickname)
% loads the detected spike waveforms (Spkwave), their spatial extent
% (Spkmask) and the noise covariance matrix (WhiteningMatAllSpatial) from
% the Spkdata file with the suffixe fnickname.
    obj.InitSpkdata();
    
    if nargin<2 || isempty(fnickname)
        fname=[obj.datafolder filesep obj.filename '_Spkdata.mat'];
    else
        fname=[obj.datafolder filesep obj.filename '_Spkdata_' fnickname '.mat'];
    end
    if exist(fname)
        load(fname,'info');
        obj.Params.detect.FcHigh=info.detect.FcHigh;
        obj.Params.detect.FcLow=info.detect.FcLow;
        obj.Params.detect.zthresh=info.detect.zthresh;            
        obj.Params.detect.wtime1=info.detect.wtime1;
        obj.Params.detect.wtime2=info.detect.wtime2;
        obj.Params.detect.SpkADgain=info.detect.SpkADgain;
        obj.Params.detect.Tcensoredfactor=info.detect.Tcensoredfactor;
        obj.Params.fileinfo.electrode=info.fileinfo.electrode;
        if isfield(info.detect,'FGroundtruth')
            obj.Params.detect.FGroundtruth = info.detect.FGroundtruth;
        else
            obj.Params.detect.FGroundtruth = false;
        end
        if isfield(info.detect,'filepathgroup')
            if numel(info.detect.filepathgroup)>1
                filename = obj.Params.detect.filepathgroup{1};
                obj.Params.detect.filepathgroup = info.detect.filepathgroup;
                obj.Params.detect.filepathgroup{1} = filename;
            end
        else
            obj.Params.detect.filepathgroup = obj.datafolder;
        end
        
        
        
        load(fname,'spkevent');
        obj.Spkevent=spkevent;
        try
            load(fname,'spkwave');
            obj.Spkwave=spkwave;
        catch
            warning('Spk waveform ensemble probably too large to be loaded');
            obj.Spkwave=[];
        end
        
        load(fname,'whiteningmatallspatial');
        obj.WhiteningMatAllSpatial=whiteningmatallspatial;
        
        load(fname,'rawdata');
        obj.RawData=rawdata;
        for ch=1:size(obj.RawData,1)
            obj.detectthresh(ch)=obj.Params.detect.zthresh*median(abs(single(obj.RawData(ch,:))/obj.Params.detect.SpkADgain)/0.6745);
        end
        load(fname,'spkmask');
        obj.Spkmask=spkmask;
    else
        warning('no spike data saved yet');
        obj.Spkwave=[];
    end
end