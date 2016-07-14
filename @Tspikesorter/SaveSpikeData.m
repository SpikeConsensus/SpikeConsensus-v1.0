function obj=SaveSpikeData(obj,fname)
% SaveSpikeData saves under the same folder as the rawdata 
% (filename_Spkdata_fname.mat) a mat file containing the detected 
% spike waveforms (Spkwave), their corresponding spatial extent 
% (Spkmask) and the noise covariance matrix (WhiteningMatAllSpatial).
    if nargin<2
        fname=[obj.datafolder filesep obj.filename '_Spkdata.mat'];
    else
        fname=[obj.datafolder filesep obj.filename '_Spkdata_' fname '.mat'];
    end
    info=obj.Params;
    spkevent=obj.Spkevent;
    spkwave=obj.Spkwave;
    spkmask=obj.Spkmask;

    whiteningmatallspatial=obj.WhiteningMatAllSpatial;
    rawdata=int16(obj.RawData(:,1:floor(0.5*end))*obj.Params.detect.SpkADgain);
    save(fname,'spkevent','spkwave','spkmask','whiteningmatallspatial','rawdata','info','-v7.3');
end