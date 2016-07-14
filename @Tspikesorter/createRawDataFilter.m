function obj=createRawDataFilter(obj)
%creates the Butterworth filters used to band-pass filter the raw data.
    Wn=[obj.Params.detect.FcLow obj.Params.detect.FcHigh]/(obj.Params.fileinfo.samplingrate/2);
    [obj.butterHB,obj.butterHA]=butter(obj.Params.detect.butterOrder,Wn(1),'high');
    [obj.butterLB,obj.butterLA]=butter(obj.Params.detect.butterOrder,Wn(2),'low');            
end